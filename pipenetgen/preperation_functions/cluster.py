import osmnx as ox
import networkx as nx
from sklearn.cluster import SpectralClustering
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import colormaps

from . import network_modification as nm
from .. import linear_programs as lp

import shapely as sp



def run_cluster_algorithm(G,n,parameters,show=False,cluster_params=None,target_topo = {'link_density':100,'node_degree':100,'meshedness':100}):
    parameters = parameters.copy()

    for edge in G.edges: 
        G.edges[edge]['roughness'] = parameters['C']

    for node in G.nodes: 
        try:
            G.nodes[node]['y'] = G.nodes[node]['lat']
            G.nodes[node]['x'] = G.nodes[node]['lon']
        except:
            G.nodes[node]['lat'] = G.nodes[node]['y']
            G.nodes[node]['lon'] = G.nodes[node]['x']
            
    subgraphs, topo, connectivity, colors, coords = make_clusters(G,n,show=show)
    print("Clusters Are Connected",check_connected(subgraphs))
    CG = prepare_cluster_graph(subgraphs,colors,coords,connectivity)

    # prepare for cluster run

    # Swap elevations, running the CG with no elevation sicne focus is soley on a theoretical flow distribution
    for node in CG.nodes():
        CG.nodes[node]['true_elevation'] = CG.nodes[node]['elevation']
        CG.nodes[node]['elevation'] = 0
    
    cg_lp = lp.milp(CG,parameters,network_type='simple')# Create entwork, simple denotes all pipes must remain
    
    seed = parameters.get('seed',0)
    cg_lp.m.params.seed = seed

        # Run
    cg_lp.load_model() 
    
    cg_lp.run_model()

    # Replace centroid elevations
    for node in CG.nodes():
        CG.nodes[node]['elevation'] = CG.nodes[node]['true_elevation']

    if show:
        cg_lp.super_plot(cg_lp.G,node_label="subgraph",edge_label="flowrate") # Showcase of flow


    run_order =  find_cluster_run_order(cg_lp.G) # Order in which centroids clusters will be run
    pseudo_demand = find_cluster_adjacency(cg_lp.G) # The theoretical flow between clusters
    print(pseudo_demand,run_order)
    # Any changes to parameters are made here
    if cluster_params:
        parameters.update(cluster_params)


    # Add in real/target values if you know them, otherwise choose overall topo,
    real_topo = target_topo

    # # Set manually

    run_sg,edge_info = run_clusters(subgraphs,CG,parameters,run_order,pseudo_demand,real_topo, seed = seed)

    base = stitch_subgraphs(G.to_undirected(),run_sg,edge_info,True)

    return base

def cluster(G,n,show=True):
    """Create n clusters in the network G

    Args:
        G (osmnx graph): The graph to cluster into smaller subgraphs
        n (int): The number of clusters
        show (bool, optional): Whether to plot the clusters. Defaults to True.

    Returns:
        list: A list of osmnx graphs (subgraphs of the orignal network)
        list: A list of labels for each node in G indicating which cluster to attach to
        list: A list of all cluster centroid locations
        list: A list of all colors used in plotting.
    """

    # Extract node coordinates from the network
    node_data = G.nodes(data=True)
    node_coords = np.array([[data['x'], data['y']] for _, data in node_data])

    # Specify the number of clusters (neighborhoods)
    n_clusters = n # Increase until all subgraphs are connected

    # Use Spectral Clustering to cluster nodes into neighborhoods
    spectral = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=0)
    labels = spectral.fit_predict(node_coords)

    # Find the centroid of each cluster
    centroids = []
    for cluster_label in range(n_clusters):
        cluster_indices = np.where(labels == cluster_label)[0]
        cluster_coords = node_coords[cluster_indices]
        centroid = np.mean(cluster_coords, axis=0)
        centroids.append(centroid)

    centroids = np.array(centroids)

    # Visualize the clustering results and centroids
    if show:
        fig, ax = ox.plot_graph(G, node_color=labels, node_size=30,  bgcolor='k',show=False)
        ax.scatter(centroids[:, 0], centroids[:, 1], c='red', marker='X', s=100, label='Centroids')
        plt.legend()
        plt.show()


    viridis_colormap = colormaps.get_cmap('viridis')

    # Get the RGBA values of n colors from the colormap
    colors_rgba = viridis_colormap(np.linspace(0, 1, n_clusters))
    # colors_rgba

    subgraphs = {}
    for cluster_label in range(n_clusters):
        # all the nodes meeting the 
        cluster_indices = np.where(labels == cluster_label)[0]
        
        # Extract subgraph based on node indices
        subgraph_nodes = G.subgraph(np.array(G.nodes)[labels==cluster_label])#G.subgraph(cluster_indices)
        
        # Add subgraph to the dictionary
        subgraphs[cluster_label] = subgraph_nodes
    # for subgraph in subgraphs.values():
    #     ox.plot_graph(subgraph)

    return subgraphs, labels, centroids, colors_rgba


def check_connected(subgraphs):
    """Check whether all subgraphs in a list are connected
    Subgraphs need to be connected to be run in the MILP

    Args:
        subgraphs (list): A list of osmnx graphs

    Returns:
        bool: Whether they are all connected
    """
    for k,s in subgraphs.items():
        if len(s.nodes) < 1:
            print("This Subgraph has no nodes")
        sg = s.to_undirected()
        if not nx.is_connected(sg):
            return False
    return True


def reassign_component(graph, component, labels, target_label):
    """Reassigns all nodes in a component list to a target label.

    Args:
        graph (osmnx graph): A graph with components
        component (list): A list of nodes that need to be reassigned 
        labels (list): A list of node cluster labels
        target_label (int): the target cluster for the components

    Returns:
        list: A new label list
    """
        
        
    for node in component:
        labels[list(graph.nodes).index(node)] = target_label
    return labels

def check_and_reconnect_subgraphs(graph, subgraphs, labels):
    """
    Checks if each subgraph is connected, finds the disconnected sections,
    and attaches the smaller of them to one of their neighboring clusters.

    Args:
        graph (networkx.Graph): The main graph.
        subgraphs (dict): A dictionary where keys are cluster labels and values are subgraphs.
        labels (list): A list of cluster labels corresponding to each node in the graph.

    Returns:
        dict: Updated subgraphs dictionary.
        list: Updated labels list.
    """

    updated_subgraphs = subgraphs.copy()
    updated_labels = labels.copy()
    
    for cluster_label, subgraph in subgraphs.items():

        if len(subgraph.nodes) == 0:  # Check if the subgraph is empty
                print(f"Warning: Cluster {cluster_label} is empty.")
                continue  # Skip empty clusters


        if not nx.is_weakly_connected(subgraph.to_directed()):
            disconnected_components = list(nx.weakly_connected_components(subgraph.to_directed()))
            
            # Find the biggest component and keep as is
            sizes = []
            for component in disconnected_components:
                sizes.append(len(component))
            
            for component in disconnected_components:
                if len(component) != max(sizes):
                    for node in component:
                        # Find the neighbors of the node in the original graph
                        neighbors = nx.neighbors(graph, node)
                        
                        if neighbors:
                            # Find the cluster of the neighbor
                            neighbor_clusters = [updated_labels[list(graph.nodes).index(neighbor)] for neighbor in neighbors]
                            
                            # Exclude the current cluster label from the neighbor clusters
                            neighbor_clusters = [nc for nc in neighbor_clusters if nc != cluster_label]
                            print(neighbor_clusters)
                            if neighbor_clusters:
                                print("Neighbor")
                                # Reassign the node to one of the neighbor's clusters
                                target_label = neighbor_clusters[0]  # Choose the first neighbor's cluster as target
                                
                                updated_labels = reassign_component(graph, component, updated_labels, target_label)
                                
                                break
                            else:
                                print("No neighbor")
                                print(node,graph.nodes[node],graph.edges(node))
    # print("Here")
    # Update subgraphs based on new labels
    updated_subgraphs = {}
    for cluster_label in set(updated_labels):
        cluster_indices = [i for i, label in enumerate(updated_labels) if label == cluster_label]
        subgraph_nodes = [list(graph.nodes)[i] for i in cluster_indices]
        subgraph = graph.subgraph(subgraph_nodes)
        updated_subgraphs[cluster_label] = subgraph
    
    return updated_subgraphs, updated_labels


def plot_clusters(G,labels):
    """ Visualize the clustering results and centroids

    Args:
        G (osmnx graph): the graph to plot
        labels (list): the cluster node labels
    """
   
    fig, ax = ox.plot_graph(G, node_color=labels, node_size=30,  bgcolor='k',show=False)
    plt.legend()
    plt.show()


def attach_cluster_labels(G,labels):
    """Attach the cluster labels to the network G 

    Args:
        G (osmnx graph): The network to receive node labels
        labels (list): The cluster labels list
    """
    for i in range(len(labels)):
        G.nodes[list(G.nodes)[i]]['cluster'] = labels[i]


def distance(p1,p2):
    """Euclidian distance

    Args:
        p1 (tuple): x,y of point 1
        p2 (tuple): x,y of point 2

    Returns:
        float: the distance
    """
    return ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**.5


def find_to_whom(G):
    """ Finds all subgraph nodes and which nodes they are connected to in other subgraphs

    Args:
        G (osmnx graph): The graph to check
    """
   
    
    for node in G.nodes: # Find all other subgraphs the nodes are connected to
        G.nodes[node]['connected'] = []
        cluster = G.nodes[node]['cluster']
        nbrs = nx.neighbors(G,node) # neighbors
        for nbr in nbrs:
            nbr_cluster =  G.nodes[nbr]['cluster'] 
            if nbr_cluster != cluster: # If it is a new cluster (border nodes)
                G.nodes[node]['connected'].append(nbr_cluster) # The cluster of my neighbor
                if "to_whom" not in G.nodes[node]:
                    G.nodes[node]['to_whom'] = []
                G.nodes[node]['to_whom'].append(nbr)

def find_connected_clusters(subgraphs,cluster_data,colors_rgba,show=True):
    """Finds all nodes that are connected to neighbors that are not in their own cluster
    Plots them
    finds the maximum topo statistics for each node

    Args:
        subgraphs (list): a list of osmnx graphs
        cluster_data (list): labels
        colors_rgba (list): all colors used for each cluster
        show (bool, optional): whether to plot. Defaults to True.

    Returns:
        dict: A dictionary of topographci statistics for each subgraph
        dict: Which clusters are connected 
        """
    connected_clusters = {i:set() for i in range(len(subgraphs))} # clusterid:[clusters_it_is_connected_to]
    sg_topo_data = {}

    for key, subgraph in subgraphs.items():
        if show:
            print("Subgraph",key)
        nc = []
        ns = []
        for n,data in subgraph.nodes(data=True):
           
            if data.get('connected'):
                 # if it has a neighbor in another cluster
                for i in range(len(data['connected'])):
                    connected_clusters[key].add(data['connected'][i])
                nc.append(colors_rgba[data['connected'][0]])
                ns.append(100)
                subgraph.nodes[n]['color'] = colors_rgba[data['connected'][0]]
            else:
                nc.append(colors_rgba[key])
                ns.append(15)
                subgraph.nodes[n]['color'] = colors_rgba[key]
            data['name'] = n
        # Max connectivity constraints
        ed = len(subgraph.edges)
        no = len(subgraph.nodes)
        
        if len(subgraph.nodes) > 1:
            link_density=2*ed/(no*(no-1))
            node_degree = 2*ed/no
            meshedness = (ed-no+1)/(2*no-5)
            sg_topo_data[key] = {"link_density":link_density,'node_degree':node_degree,'meshedness':meshedness}
            ox.plot_graph(subgraph,node_color=nc,node_size=ns,show=show,close=not show)
    return sg_topo_data,connected_clusters

def prepare_cluster_graph(subgraphs,colors_rgba,centroids,connected_clusters,show=True):
    """Prepare to run the MILP on all centroids. Where the centroid represents the total demand of the network.

    Args:
        subgraphs (list): A list of osmnx subgraphs
        colors_rgba (list): A list of colors of subgraphs
        centroids (list): Centroid locations
        connected_clusters (dict): Which subgraphs are connected and how much
        show (bool, optional): Whether to plot. Defaults to True.

    Returns:
        osmnx graph: A graph ready to optimize.
    """
    CG = nx.MultiDiGraph()

    # Add nodes with positions to the graph
    for i in range(len(centroids)):
        if len(subgraphs[i]) > 0: # Not blank
            elev = np.mean(list(nx.get_node_attributes(subgraphs[i],'elevation').values()))
            dem = sum(list(nx.get_node_attributes(subgraphs[i],'demand').values()))
            CG.add_node(i, pos=(centroids[i][0],centroids[i][1]),
                        x = centroids[i][0],
                        y=centroids[i][1],
                        color=colors_rgba[i],
                        subgraph=i,
                        elevation=elev,
                        demand = dem)

    for key,connections in connected_clusters.items():
        for connection in connections:
            if key and connection in CG.nodes:
                length = distance(CG.nodes[key]['pos'],CG.nodes[connection]['pos'])
                CG.add_edge(key,connection,length=length,exists='both')

    CG.graph['crs'] = subgraphs[0].graph['crs']

    # Plot the graph
    pos = nx.get_node_attributes(CG, 'pos')
    n_demands = nx.get_node_attributes(CG,"demand")
 
    nx.draw(CG, pos, with_labels=True, node_size=700, node_color=colors_rgba, font_size=8, font_color='black', font_weight='bold')
    
    # Display the plot
    if show:
        plt.show()
    else:
        plt.close()
    return CG


def find_cluster_run_order(G):
    """From the output of the Cluster MILP find an order that ensures no cluster is run before its predecessors has been run.

    Args:
        G (osmnx graph): a directed osmnx graph showcasing how water flows through the network

    Returns:
        list: A list of what order to run each subgraph
    """
    run_order = list(nx.topological_sort(G))
    return run_order

def find_cluster_adjacency(G):
    """Make a flow adjacency matrix (dict) for a clustered graph 

    Args:
        G (osmnx graph): Directed graph showing flow direction between centroids

    Returns:
        dict(dict): A dictionary of nodes with values of dictionaries of keys neighbor nodes and their flow as values  
    """
    # assign nodal pseudo demands in an adjacency matrix
    # it will be the same as flow in since this is enough to meet its demands and neighbors
    pseudo_demand = {}
    for edge in G.edges:
        u,v,k = edge
        if u not in pseudo_demand:
            pseudo_demand[u] = {}
        if v not in pseudo_demand:
            pseudo_demand[v] = {}
        if v not in pseudo_demand[u]:
            pseudo_demand[u][v] = 0
        if u not in pseudo_demand[v]:
            pseudo_demand[v][u] =0

        pseudo_demand[u][v] = round(G.edges[edge]['flowrate'],6)
        pseudo_demand[v][u] = round(-1 * G.edges[edge]['flowrate'],6)


    return pseudo_demand


def make_clusters(G,n_cluster,show=True): 
    """Run all steps of the clusterig algorithm at once to produce a G that has been split up"""
    # Group into clusters, they may not be hydraulicly connected so reconnect
    subgraphs, labels, centroids, colors_rgba = cluster(G,n_cluster,show) 

    # reconect to ensure all nodes in clusters are connected
    subgraphs,labels = check_and_reconnect_subgraphs(G.to_undirected(),subgraphs, labels)


    if show:
        plot_clusters(G,labels) # Now those areas have been reconnected
        plt.show()
    attach_cluster_labels(G,labels) # Correctly reattach labels

    # Attach data back to G
    G = G.to_undirected()
    cluster_data = nx.get_node_attributes(G,'cluster')

    find_to_whom(G) # Find which edges are attached to other edges

    subgraphs = {} # Split up G into subgraphs by their label
    for i in range(n_cluster):
        subgraphs[i] = G.subgraph(np.array(G.nodes)[labels==i])

    sg_topo_data,connected_clusters = find_connected_clusters(subgraphs,cluster_data,colors_rgba,show=show) # Gets nodes with other clusters and then attaches that info to the node
    return subgraphs, sg_topo_data,connected_clusters,colors_rgba,centroids

def run_clusters(subgraphs,CG,parameters,run_order,pseudo_demand,real_topo,seed = 0):
    # Does LP iteration of subgraphs
    run_sg = {}

    big_nodes = [] # When a node has a big pipe into it add it to this list, Use this when checking clusters connections to see if should group
    clusters_run = [] # When a cluster has been run
    edge_info = {}
    sg_topo_data = {k:{} for k in subgraphs.keys()}


    # for sgoi,sg in subgraphs.items():
    for sgoi in run_order:

        # Prep subgraph
        sg = subgraphs[sgoi]
        sg = sg.copy()
        nm.rename_edge_keys(sg)
        net = sg.copy()
        net = net.to_directed()

        # Assign flow from cluster
        for cluster in pseudo_demand[sgoi]:
            net.add_node("C"+str(cluster),**CG.nodes[cluster])
            net.nodes["C"+str(cluster)]['demand'] = -1*pseudo_demand[sgoi][cluster]


        eoi = {}
        noi = []
        # Ensure that clusters are connected by a large pipe using a group and anti group method
        # Also ensures that pipes don't pull large pipes out of nowhere. 
        # Loosely coupled clusters
        for node in net.nodes:
            if "to_whom" in net.nodes[node]:
                tw,c = net.nodes[node]['to_whom'],net.nodes[node]['connected']
                noi.append(node)
                for i in range(len(tw)):
                    length = 10 

                    grp = None
                    anti_grp = None
                    if pseudo_demand[sgoi][c[i]] <= 0: # Going to a demand cluster
                        grp = c[i]
                        eoi[(node,"C"+str(c[i]))] = (node,tw[i])
                    else:
                        eoi[("C"+str(c[i]),node)] = (tw[i],node)
                    if tw[i] not in big_nodes:
                        anti_grp = c[i]
                    
                    net.add_edge(node,"C"+str(c[i]),0,length=length,og=(node,tw[i]),group = grp)
                    net.add_edge("C"+str(c[i]),node,0,length=length,og=(tw[i],node),anti_group = anti_grp)

        # Calculate Topo Constraints
        ed = len(net.edges)/2 
        no = len(net.nodes)
        link_density=2*ed/(no*(no-1))
        node_degree = 2*ed/no
        meshedness = (ed-no+1)/(2*no-5)
        
        sg_topo_data[sgoi]['link_density'] = link_density
        sg_topo_data[sgoi]['node_degree'] = node_degree
        sg_topo_data[sgoi]['meshedness'] = meshedness
        print(sg_topo_data[sgoi])
        print(real_topo)
        
        # Ensure topo input isn't bigger than possible
        for key,val in sg_topo_data[sgoi].items():
            if val >= real_topo[key]: # if the maximum possible amount is bigger than the desired amount use desired amount
                val = real_topo[key]
            else: # If the real amount is above the maximum use the max amount
                val = val
            parameters[key] = int(val *10000)/10000 # Round down and then multiply
        print(parameters['link_density'],parameters['node_degree'],parameters['meshedness'])
        # Plot original
        ox.plot_graph(net)
    
        nm.node_elev_mult(net,parameters['elev_mult']) # Correct elevations

        # Run
        LP_model2 = lp.milp(net,parameters,network_type="cluster")
        LP_model2.load_model()
        LP_model2.m.params.seed = seed
        LP_model2.run_model()
    
        # Save to a dict and plot
        run_sg[sgoi] = LP_model2.G
        LP_model2.plot_network(LP_model2.G,"diameter","flowrate")

        # Save info about which nodes have big ability
        for node in noi:
            for edge in LP_model2.G.edges(node):
                u,v = edge
                diam_in = LP_model2.edge_diam[(u,v)][0.6].x
                diam_out = LP_model2.edge_diam[(v,u)][0.6].x
                if diam_in+diam_out >= 1:
                    big_nodes.append(node)

        for edge,og in eoi.items(): # Save information about edges diameters and flows
            if edge in LP_model2.G.edges:
                u,v = edge
                diam = LP_model2.actual_diam[(u,v)].x
                flow = LP_model2.edge_flow[(u,v)].x
                edge_info[(og[0],og[1])] = {"diameter":diam,"flowrate":flow}

    return run_sg,edge_info


def calculate_geometry(G,node1,node2):
    """ Calcualtes straight geometery"""
    p1 = sp.Point(G.nodes[node1]['x'],G.nodes[node1]['y'])
    p2 = sp.Point(G.nodes[node2]['x'],G.nodes[node2]['y'])
    return sp.LineString([p1,p2])


def stitch_subgraphs(G,run_sg,edge_info,show=True): 
    base = run_sg[0].copy()

    for graph in range(1,len(run_sg)):
        base = nx.compose(base,run_sg[graph])

    base_raw = base.copy()


    nm.rename_edge_keys(G) # For adding geoms need keys to be 0

    eoi = []
    for edge,data in edge_info.items():
        u,v = edge
        u_c = G.nodes[u]['cluster']
        v_c = G.nodes[v]['cluster']
        if edge in G.edges:
            print(edge)
            data.update(G.edges[(edge[0],edge[1],0)])
        if ("C"+str(u_c),v) in base_raw.edges("C"+str(u_c)):
            base.add_edge(u,v,0,**data)
            base.edges[(u,v,0)]['geometry'] = G.edges[(u,v,0)].get('geometry',calculate_geometry(base,u,v))
            eoi.append((u,v,0))
        elif ("C"+str(v_c),u) in base_raw.edges("C"+str(v_c)):
            base.add_edge(v,u,0,**data)
            base.edges[(v,u,0)]['geometry'] = G.edges[(v,u,0)].get('geometry',calculate_geometry(base,v,u))
            eoi.append((v,u,0))


    ec = ['r' if edge in eoi else "gray" for edge in base.edges]
    ox.plot_graph(base,node_size=2,edge_color=ec,show=show,close=not show)

    base.remove_nodes_from(["C"+str(sg) for sg in run_sg.keys()] )

    ec = ['r' if edge in eoi else "gray" for edge in base.edges]
    ox.plot_graph(base,node_size=2,edge_color=ec,show=show,close=not show)

    return base


def set_reservoir_pressure_to_pump(G):
    for node,types in nx.get_node_attributes(G,"special").items():
        if types == "Reservoir":
            max_nbr = 0
            for nbr in list(nx.neighbors(G.to_undirected(),node)):
                max_head = max(abs(G.nodes[node]['head'] - G.nodes[node]['elevation']),0)# abs(G.nodes[nbr]['head'] - G.nodes[nbr]['elevation']))
                max_nbr = max(max_nbr,max_head)
            for nbr in list(nx.neighbors(G.to_undirected(),node)):
                eoi = (node,nbr,0)
                if G.nodes[nbr].get('special',False) == False: # Don't add pump between sources
                    if eoi in G.edges:             
                    
                        G.edges[eoi]['pumphead'] = max_nbr#max(abs(G.nodes[node]['head'] - G.nodes[node]['elevation']), abs(G.nodes[nbr]['head'] - G.nodes[nbr]['elevation'])) # Either the pressure, or neighbors head
                        G.edges[eoi]['pump'] = True
                        G.edges[eoi]['special'] = "pump"
                    elif (nbr,node,0) in G.edges: # If edge facing wrong way
                        data = G.edges[(nbr,node,0)] # Collect data from the edge that does exist
                        G.add_edge(eoi[0],eoi[1],0,**data) # Add a new edge out from the source
                        G.remove_edge(eoi[1],eoi[0],0) # Remove the reverse the edge
                        G.edges[(eoi[0],eoi[1],0)]['pumphead'] = max_nbr#max(abs(G.nodes[node]['head'] - G.nodes[node]['elevation']), abs(G.nodes[nbr]['head'] - G.nodes[nbr]['elevation'])) # Either the pressure, or neighbors head
                        G.edges[eoi[0],eoi[1],0]['pump'] = True
                        G.edges[eoi[0],eoi[1],0]['special'] = "pump"
            G.nodes[node]['pressure'] = 0
