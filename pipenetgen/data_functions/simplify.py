import networkx as nx
import osmnx as ox
import numpy as np
import shapely as sp

import pandas as pd

# My methods
from . import intersecting_method as im
from .. import preperation_functions as pf



def remove_road_by_type(G,types):
    """
    Removes edges of type types, UNLESS they are the only thing keeping subsections connected
    """
    G = im.trunkator(G,types)
    return G

def remove_excess_edges(G):
    if not G.is_directed():
        G = G.to_directed()
    G = im.fix_edges(G)
    return G

def search_key_dict(key_dict,value):
    out = {}
    for k,d in key_dict.items():
        out[k] =  d[value]
    return out

def consolidate_intersections(G,tol=25): 
    crs = G.graph['crs'] # Original projection
    data = {node:G.nodes[node] for node in nx.get_node_attributes(G,"special").keys()}
    G_proj = ox.project_graph(G)
    G_con = ox.consolidate_intersections(G_proj, rebuild_graph=True, tolerance=tol, dead_ends=True)

    # ADd back in sources if they got consolidated
    G_fin = ox.project_graph(G_con,to_crs = crs) # Project back
    source_fin = {node:G_fin.nodes[node] for node in nx.get_node_attributes(G_fin,"special").keys()}
    searched_fin = search_key_dict(source_fin,'name')
   
    searched_data = search_key_dict(data,'name')
  

    added_sources = []
    for k,d in data.items():
        if d['name'] not in list(searched_fin.values()):
            print(d['special'])
            G_fin.add_node(k,**d)
            print("Added",d['name'],'as',k)
            added_sources.append(k)
            acceptable_connections = list(set(G_fin.nodes) - set(source_fin.keys()) - set(added_sources)) # This list isn't working. It still connects new sources to existing sources
            # print(set(source_fin.keys()),acceptable_connections)
            pf.network_modification.connect_to_closest(G_fin,k,acceptable_connections)


    return G_fin

def remove_excess_nodes(G,iterations=25):
    G,rem = im.fix_nodes(G)
    while rem != 0 and iterations >= 0:
        G,rem = im.fix_nodes(G)
        iterations-=1
    return G

def remove_multi_edges(G,show=False):
    """Remove Multi-edges"""
    seen = set()
    edges_to_remove = []

    # G = nx.convert_node_labels_to_integers(G)

    for edge in G.edges:
        u,v,k = edge
        tup = max(u,v),min(u,v)
        if tup in seen:
            edges_to_remove.append(edge)
        else:
            seen.add(tup)

    print(len(edges_to_remove))
    G.remove_edges_from(edges_to_remove)
    ec = ['r' if edge in edges_to_remove else "y" for edge in G.edges]
    ox.plot_graph(G,node_size=2,edge_color=ec,show=show,close=not show)
    return G

def fix_geometry(G,show): # TODO: Look into this for the network to inp function (also see whether wntr can do that too)
    #If geometry is execsively large (like doubles up on itself, replace with straight)
    fixed_edges = []

    mult = 0.5 #This isn't perfect

    for edge in G.edges:
        u,v,k = edge
        
        if 'geometry' in G.edges[edge]:
            length = G.edges[edge]['geometry'].length
        else:
            G.edges[edge]['geometry'] = im.calculate_geometry(G,edge[0],edge[1])
            length = G.edges[edge]['geometry'].length
        
        coords2 = list(G.edges[edge]['geometry'].coords)[1:-1]

        if len(coords2) > 2: #Captures the issue where consolidate intersections rebuilds badly and flips sides
            n1,n2 = im.calculate_geometry(G,edge[0],edge[1]).coords
            geo2 = sp.LineString([n1]+coords2[::-1]+[n2])
            length2 = geo2.length

        straight = sp.LineString([G.edges[edge]['geometry'].coords[0],G.edges[edge]['geometry'].coords[-1]]).length
        dif = length-straight

        if dif > mult*length:
            if length2-straight<dif:
                G.edges[edge]['geometry'] = geo2
                G.edges[edge]['length'] = length2
            else:
                G.edges[edge]['geometry'] = im.calculate_geometry(G,edge[0],edge[1])
                G.edges[edge]['length'] = G.edges[edge]['geometry'].length
            
            fixed_edges.append(edge)

    ec = ['r' if edge in fixed_edges else "y" for edge in G.edges]
    fig,ax = ox.plot_graph(G,node_size=1,edge_linewidth=1,edge_color=ec,show=show,close=not show)
    return G


def save_og_name(G):
    for node in G.nodes:
        name = G.nodes[node].get("Name",G.nodes[node].get("name",str(node)))
        G.nodes[node]['name'] = name

def fix_topo(G,show=False,consolidate_dist=25,types=['trunk','motorway','motorway_link','trunk_link','tertiary_link','primary_link']):
    if G.graph.get("topo_fixed"):
        print("Topography has already been fixed, cannot re-consololidate instersections")
    else:
        save_og_name(G)
        G = nx.convert_node_labels_to_integers(G)
        G = remove_road_by_type(G,types)
        G = remove_excess_edges(G)
        G = consolidate_intersections(G,tol=consolidate_dist)
        G = remove_excess_nodes(G)
        G = remove_excess_edges(G)
        G = remove_multi_edges(G,show)
        G = fix_geometry(G,show)
        G.graph['topo_fixed'] = True
    return G

    
    