# Purpose: To receive an osmnx network and return one where there are no intersecting edges
import networkx as nx
import matplotlib.pyplot as plt
import shapely as sp
import geopandas as gp
import numpy as np
import osmnx as ox
import pandas as pd
# from sklearn.cluster import DBSCAN
# from . import network_to_inp as ntoi



def trunkator(G,types):
    """Removes edges in types unless they are needed to connect graph"""
    G2 = G.copy()
    prim = [] #Get all primary roads 
    for edge in G2.edges:
        if G2.edges[edge].get('highway') in types:
            prim.append(edge)
    ec = ["y" if edge in prim else "r" for edge in G2.edges] 
    # fig,ax = ox.plot_graph(G2,edge_color=ec,node_size = 2)
    # plt.show()

    G2 = G2.to_undirected()
    Gpt = G2.copy()
    copied_graph = G2.copy()
    # Remove the specified edges from the copied graph
    copied_graph.remove_edges_from(prim)

    # Find connected components in the copied graph
    connected_components = list(nx.connected_components(copied_graph))
    # h_node = max(connected_components)[0] # a node in the largest component group
    h_node= list(max(connected_components))[0]


    n_to_e = {}
    for edge in prim:
        u,v,k = edge
        for n in u,v:
            if n in n_to_e:
                n_to_e[n].append(edge)
            else:
                n_to_e[n] = [edge]


    edges_to_keep = []

    # Iterate through the connected components and check for connectivity
    for component in connected_components:
        if component != max(connected_components):
            length = len(component)
            if length > 0: #Decide that only if it is big then is it important
                target = list(component)[0]
                
                paths = ox.routing.k_shortest_paths(Gpt, h_node, target, 1, weight='length')
                
                for path in paths:
                    for node in range(len(path)-1):
                        n = path[node]
                        if n in n_to_e:
                            edges = n_to_e[n]
                            edges_to_keep.extend(edges)

    #print(len(prim),len(edges_to_keep))

    edges_to_remove = list(set(prim) - set(edges_to_keep))
    print("removing",len(edges_to_remove))

    ec = ['y' if edge in edges_to_keep else 'r' if edge in edges_to_remove else "w" for edge in Gpt.edges]

    G4 = G2.copy()
    G4.remove_edges_from(edges_to_remove)
    return G4


def calculate_geometry(G,node1,node2):
    """ Calcualtes straight geometery"""
    p1 = sp.Point(G.nodes[node1]['x'],G.nodes[node1]['y'])
    p2 = sp.Point(G.nodes[node2]['x'],G.nodes[node2]['y'])
    return sp.LineString([p1,p2])


def fix_nodes(G):
    """If they have no neighbors get rid of em"""
    G_copy = G.copy()
    #look at count_streets_per_node if it is <= 1 then node is removed.
    nodes_to_remove = []
    for node,val in ox.stats.count_streets_per_node(G_copy).items():
        if val <= 1:
            if not G.nodes[node].get('special'):
                nodes_to_remove.append(node) 
    print("Removing",len(nodes_to_remove),"nodes")
    G_copy.remove_nodes_from(nodes_to_remove)
    return G_copy,len(nodes_to_remove)


def fix_edges(G):
    """
    1. Will add geometry and length to edges
    2. Will remove edges that:
     - Are duplicates from point A to B (keeps shortest path)
     - go from point A to A (loops)
    returns G2
    """

    edges_to_remove = []

    prexisting = {}
    G_copy = G.copy()

    for edge in G_copy.edges:
        u,v,n = edge
        
        if 'geometry' not in G_copy.edges[edge]: #Add geometry
            geometry = calculate_geometry(G_copy,u,v)
            G_copy.edges[edge]['geometry'] = geometry
        if 'length' not in G_copy.edges[edge]: #Add length
            edge_length = geometry.length
            G_copy.edges[edge]['length'] = edge_length
        edge_length = G_copy.edges[edge]['length']

        if type(u) != type(v):
            u = str(u)
            v = str(v)
        tupl = (max(u,v),min(u,v))#Create a non unique identifier that connects nodes
        

        if u == v: #If it self loops
            edges_to_remove.append(edge)
        else:
            if tupl not in prexisting: #See if we have logged it before
                prexisting[tupl] = (edge_length,edge)
            
            elif edge_length <= prexisting[tupl][0]: #if it has shorter length then remove the current edge and log this as new edge
                edges_to_remove.append(prexisting[tupl][1])
                prexisting[tupl] = (edge_length,edge)
    
    print("Removing", len(edges_to_remove), "excess Edges")

    G_copy.remove_edges_from(edges_to_remove)
    return G_copy
