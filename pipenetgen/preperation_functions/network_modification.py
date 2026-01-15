
import networkx as nx
import osmnx as ox
import pandas as pd
import numpy as np
import requests
import math

def assign_proportional_demand(G):
    td = get_total_node_demands(G)
    for node in nx.get_node_attributes(G,"special"):
        td -= G.nodes[node]['demand']
    for node in nx.get_node_attributes(G,"special"):
        G.nodes[node]['demand'] = -G.nodes[node]['Proportion'] * td


def assign_meteo_elevation(G,crs="EPSG:4326",grp_size=5):
    """Using meteo elevation api get elevation data for all nodes. Less precise than tif files. Less storage though.
    If your graph is already projected OR the x,y coordinates are already lon,lat. Then use None for crs

    Args:
        G (osmnx graph): the graph to apply elevations to
        crs (str, optional): the CRS to apply your graph to. Defaults to 'EPSG:4326'.
    """
    # Project graph
    if crs:
        G_proj = ox.project_graph(G,to_crs=crs)

    # Capture data 
    data = {}
    node_set = {}
    n = 0
    i = 1

    no_lat_nodes = [] # Don't have geo info (shouldn't exist)

    # Every 10 nodes make a group
    for node in list(G_proj.nodes):
        if n not in node_set:
            node_set[n] = []
        if 'x' in G_proj.nodes[node]:
            node_set[n].append(node)
        else:
            print("missing")
            no_lat_nodes.append(node)
        if i == grp_size:
            i = 0
            n += 1
        i += 1
    
    # Make the request based on node groups
    for grp,nodes in node_set.items():
        lats = []
        lons = []
        for node in nodes:
            if 'x' in G_proj.nodes[node]:
                lat = G_proj.nodes[node]['y']
                lon = G_proj.nodes[node]['x']
                lats.append(str(lat))
                lons.append(str(lon))
                data[node] = {'y':lat,"x":lon}
            else:
                print("HEY")
    
        response = requests.get(f"https://api.open-meteo.com/v1/elevation?latitude={','.join(lats)}&longitude={','.join(lons)}")
        # print(response)
        response_list = response.json()['elevation']
        # print(response.json())
        if response.status_code == 200:
            for i in range(len(nodes)):
                node = nodes[i]
                # print(response_list[i]['latitude'] , G_proj.nodes[node]['y'])
                G.nodes[node]['elevation'] = response_list[i]
        time.sleep(0.01) # To avoid rate limiting

import time
def assign_open_elevation(G,crs='EPSG:4326'):
    """Using open elevation api get elevation data for all nodes. Less precise than tif files. Less storage though.
    If your graph is already projected OR the x,y coordinates are already lon,lat. Then use None for crs

    Args:
        G (osmnx graph): the graph to apply elevations to
        crs (str, optional): the CRS to apply your graph to. Defaults to 'EPSG:4326'.
    """
    # Project graph
    if crs:
        G_proj = ox.project_graph(G,to_crs=crs)

    # Capture data 
    data = {}
    node_set = {}
    n = 0
    i = 1

    no_lat_nodes = [] # Don't have geo info (shouldn't exist)

    # Every 10 nodes make a group
    for node in list(G_proj.nodes):
        if n not in node_set:
            node_set[n] = []
        if 'x' in G_proj.nodes[node]:
            node_set[n].append(node)
        else:
            print("missing")
            no_lat_nodes.append(node)
        if i == 4:
            i = 0
            n += 1
        i += 1
    
    # Make the request based on node groups
    for grp,nodes in node_set.items():
        locations = []
        for node in nodes:
            if 'x' in G_proj.nodes[node]:
                lat = G_proj.nodes[node]['y']
                lon = G_proj.nodes[node]['x']
                locations.append(f"{lat},{lon}")
                data[node] = {'y':lat,"x":lon}
            else:
                print("HEY")
        
        response = requests.get(f"https://api.open-elevation.com/api/v1/lookup?locations={'|'.join(locations)}")
        response_list = response.json()['results']
        # print(response.json())
        time.sleep(0.01) # To avoid rate limiting
        if response.status_code == 200:
            for i in range(len(nodes)):
                node = nodes[i]
                # print(response_list[i]['latitude'] , G_proj.nodes[node]['y'])
                G.nodes[node]['elevation'] = response_list[i]['elevation']


def find_lowest_node(G):
    source, elev =  min(list(nx.get_node_attributes(G,'elevation').items()),key=lambda x: x[1])
    return (source,elev)

def find_furthest_node(graph, allowable_nodes, selected_nodes):
    furthest_node = None
    max_distance = -1
    
    for node in allowable_nodes:
        if node in selected_nodes:
            continue
        # Calculate the minimum distance from this node to any of the selected nodes
        min_distance = min(nx.shortest_path_length(graph, source=node, target=sn, weight='weight') for sn in selected_nodes)
        
        if min_distance > max_distance:
            max_distance = min_distance
            furthest_node = node
    
    return furthest_node

def find_furthest_apart_nodes(graph, n,pcent=0.4):
    # Step 1: Identify the node with the lowest elevation
    lowest_node = find_lowest_node(graph)[0]
    allowable_nodes = find_lowest_nodes(graph,nx.get_node_attributes(graph,"elevation"),pcent)
     # Initialize the list of selected nodes with the lowest node
    selected_nodes = [lowest_node]
    
    # Step 2: Iteratively find the furthest node from the selected nodes
    while len(selected_nodes) < n:
        print(selected_nodes)
        furthest_node = find_furthest_node(graph, allowable_nodes,selected_nodes)
        selected_nodes.append(furthest_node)
    
    return selected_nodes


def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def find_closest_node(G, noi, source_list):
    try:
        x_given, y_given = G.nodes[noi]['lat'], G.nodes[noi]['lon']
        mult = 111139 # Multiply geographic lat lon by this to get meters
    except:
        
        x_given, y_given = G.nodes[noi]['x'], G.nodes[noi]['y']
        mult = 1 # SOME BAD CODING #TODO: GO THROUGH ALL MENTIONS OF LAT LONG AND MAKE RIGH
    min_distance = float('inf')
    closest_node = None

    for node in G.nodes():
        try:
            x_node, y_node = G.nodes[node]['x'], G.nodes[node]['y']
            distance = euclidean_distance(x_given, y_given, x_node, y_node)
            
            if distance > 10: # I got coords wrong #TODO: THis is a terrible fix. Again with the lat lon problems.
                x_node, y_node = G.nodes[node]['y'], G.nodes[node]['x']
                distance = euclidean_distance(x_given, y_given, x_node, y_node)
            if distance < min_distance and node != noi:
                min_distance = distance
                closest_node = node
        except:
            print("Connection Error",node,noi)
    if closest_node == None:
        print("Critical Error")
    
    return closest_node,min_distance*mult # Multiply geographic lat lon by this to get meters



def connect_to_closest(G,node, source_list):
    noi,l = find_closest_node(G,node,source_list)
    # geom = png.calculate_geometry(G,node,noi)
    if G.nodes[node]['special'] == "Tank":
        G.add_edge(node,noi,length=l,highway='tank_connection')#,geometry=geom)
    else:
        G.add_edge(node,noi,length=l,highway='source_connection',special='pump',pump=True,pumphead=1)#,geometry=geom)
    # G.add_edge(noi,node,length=l)


def add_sources_from_csv(G,filename):
    df = pd.read_csv(filename)
    source_list = list(df['Name'])
    for row in df.iterrows():
        data = row[1]
        data['special'] = "Reservoir"
        data['x'],data['y'] = data['lon'],data['lat'] #TODO Fix this!!!
        # print(data["Name"])
        G.add_node(data['Name'],**data)

        connect_to_closest(G,data['Name'],source_list)

def add_tanks_from_csv(G,filename):
    df = pd.read_csv(filename)
    source_list = list(df['Name']) # TODO: Remake this nomenclature so it is lowercase name
    for row in df.iterrows():
        data = row[1]
        data['special'] = "Tank"
        data['x'],data['y'] = data['lon'],data['lat']
        print(data["Name"])
        if 'vol_curve_name' in data:
            del data['vol_curve_name']
        G.add_node(data['Name'],**data)
        
        
        connect_to_closest(G,data['Name'],source_list)

def add_sources(G,locations=None,n=1,percentage=0.4):
    """
    add sources in lowest positions (I have this coded somewhere, find it)
    """
    
    
    if locations:
        num = len(locations)
        for node in locations:
            if G.nodes[node].get("special",False): # If it already has tag, keep it (preserves tanks)
                pass
            else:
                G.nodes[node]['special'] = "Reservoir" # otherwise give it a reservoir tag
        try:
            total_dem = get_total_demands(G) - get_total_source_demands(G)
        except:
            total_dem = 1
        each_dem = total_dem/num
        for node in locations:
            G.nodes[node]['demand'] = -1* each_dem

    else:
        if n == 1:
            source, elev =  min(list(nx.get_node_attributes(G,'elevation').items()),key=lambda x: x[1])
            G.nodes[source]['demand'] = G.nodes[source]['demand'] - sum(list(nx.get_node_attributes(G,"demand").values()))
            G.nodes[source]['special'] = 'Reservoir'
        else:
            for node in G.nodes:
                G.nodes[node]['pos'] = (G.nodes[node]['x'],G.nodes[node]['y'])
            elevation_dict = nx.get_node_attributes(G,"elevation")
            source1 = min(elevation_dict,key=elevation_dict.get)


            sources = find_furthest_apart_nodes(G.to_undirected(),n,percentage)
            add_sources(G,locations=sources)

def rename_edge_keys(G):
    """Rename edge keys that are 1 when they should be 0"""
    for edge in list(G.edges):
        u,v,k = edge
        if k >= 1 and (u,v,0) not in list(G.edges):
            attrs = G.edges[edge]
            G.add_edge(u,v,0,**attrs)
            G.remove_edge(u,v,k)

def add_demands(G,demand):
    for node in G.nodes:
        G.nodes[node]['demand'] = demand


# Demand functions
def get_total_source_demands(G):
    sources = list(nx.get_node_attributes(G,"special").keys())
    total = 0
    for source in sources:
        total += G.nodes[source]['demand']
    return total

def get_total_demands(G):
    return sum(list(nx.get_node_attributes(G,"demand").values()))

def get_total_node_demands(G):
    total = 0
    for node in G.nodes:
        if G.nodes[node]['demand'] >0:
            total+=G.nodes[node]['demand']
    return total

def redistribute_demands(G):
    # Redistribute demand
    total_demand = get_total_demands(G) - get_total_source_demands(G)
    if len(nx.get_node_attributes(G,'Proportion')) == len(nx.get_node_attributes(G,'special')): # If known distrivution
        for node in nx.get_node_attributes(G,'Proportion'):
            G.nodes[node]['demand'] = -1* G.nodes[node]['Proportion'] * total_demand
        G.nodes[node]['demand'] -= get_total_demands(G)
    else: # Distribute equally
        number_of_sources = len(nx.get_node_attributes(G,'special'))
        demand = total_demand/number_of_sources

        for node in nx.get_node_attributes(G,'special'):
            G.nodes[node]['demand'] = -1 * demand
    # print(get_total_demands(G) - get_total_source_demands(G))

def node_elev_mult(G,elev_mult):
    for node in G.nodes:
        if "elevation" in G.nodes[node]:
            G.nodes[node]['elevation'] *= elev_mult



def get_minimum_positive_demand(G):
    darr = np.array(list(nx.get_node_attributes(G,"demand").values()))
    return min(darr[darr>=0])
                
def get_minimum_negative_demand(G):
    return(min(list(nx.get_node_attributes(G,"demand").values())))

def get_topo_information(G):
    ed = len(G.edges) # TODO: Should this be net?
    no = len(G.nodes)
    link_density=2*ed/(no*(no-1))
    node_degree = 2*ed/no
    meshedness = (ed-no+1)/(2*no-5)
    
    return {'link_density':link_density,"node_degree":node_degree,"meshedness":meshedness}


def correct_attr_type(G):
    for node, data in G.nodes(data=True):
        for attr, value in data.items():
            if isinstance(value, str):
                try:
                    data[attr] = int(value)
                    
                except ValueError:
                    try:
                        data[attr] = float(value)
                    except:
                        pass
                    # If conversion fails, leave the value as is
                    pass

    # Check and correct edge attributes
    for u,v, data in G.edges(data=True):
        for attr, value in data.items():
            if isinstance(value, str):
                try:
                    data[attr] = float(value)
                except ValueError:
                    # If conversion fails, leave the value as is
                    pass


def find_lowest_nodes(G, elevation_dict, percentage):
    # Sort nodes based on elevation
    sorted_nodes = sorted(elevation_dict.keys(), key=lambda x: elevation_dict[x])

    # Find the number of nodes to select based on the percentage
    num_nodes = int(len(sorted_nodes) * percentage)

    # Select the lowest nodes
    lowest_nodes = sorted_nodes[:num_nodes]

    return lowest_nodes


def give_node_pos(G,node):
    G.nodes[node]['pos'] = (G.nodes[node]['x'],G.nodes[node]['y'])

def calc_distance(G, node1, node2):
    give_node_pos(G,node1)
    give_node_pos(G,node2)
    return np.linalg.norm(np.array(G.nodes[node1]['pos']) - np.array(G.nodes[node2]['pos']))

def give_all_edges_length(G):
    for edge in G.edges:
        u,v,k = edge
        if not G.edges[edge].get('length'):
            G.edges[edge]['length'] = calc_distance(G,u,v)

def calc_all_distances(G, node, lowest_nodes):
    dists = {}
    for ln in lowest_nodes:
        dists[ln] = calc_distance(G, node, ln)

    return dists

