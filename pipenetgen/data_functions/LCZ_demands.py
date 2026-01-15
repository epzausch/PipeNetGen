import networkx as nx
import osmnx as ox
import rasterio
import geopandas as gpd
from rasterio.plot import show as raster_show

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from scipy.spatial import KDTree

def LCZ_colormap():
    data = [
        (1, 140, 0, 0, 255),
        (2, 209, 0, 0, 255),
        (3, 255, 0, 0, 255),
        (4, 191, 77, 0, 255),
        (5, 255, 102, 0, 255),
        (6, 255, 153, 85, 255),
        (7, 250, 238, 5, 255),
        (8, 188, 188, 188, 255),
        (9, 255, 204, 170, 255),
        (10, 85, 85, 85, 255),
        (11, 0, 106, 0, 255),
        (12, 0, 170, 0, 255),
        (13, 100, 133, 37, 255),
        (14, 185, 219, 121, 255),
        (15, 0, 0, 0, 255),
        (16, 251, 247, 174, 255),
        (17, 106, 106, 205, 255)
    ]


    # Extract RGBA values and normalize them to [0, 1]
    colors = [(r/255, g/255, b/255, a/255) for _, r, g, b, a in data]

    # Create a list of tuples (position, color)
    color_list = [(pos, col) for pos, col in zip(list(range(17)), colors)]

    # Create the colormap
    cmap = ListedColormap([col for _, col in color_list],N=17)
    return cmap


def load_lcz_raster(filename):
    with rasterio.open(filename) as src:
        raster_data = src.read(1)

        transform = src.transform

        crs = src.crs
    return raster_data,transform,crs


def plot_raster(raster_data,transform,**kwargs):
    raster_show(raster_data,transform=transform,vmin=1,vmax=17,**kwargs)

def make_LCZ_coords(raster_data,transform):
    rows,cols = np.indices(raster_data.shape)
    xs,ys = rasterio.transform.xy(transform,rows,cols)

    xs = np.array(xs).flatten()
    ys = np.array(ys).flatten()

    raster_coords = np.vstack((xs,ys)).T

    return raster_coords


def find_closest_node_to_LCZ(G,raster_data,transform,max_dist=0.001,show=True,ax=None):
    # Add position data to G
    for node in G.nodes:
        G.nodes[node]['pos'] = (G.nodes[node]['x'],G.nodes[node]['y'])

    node_positions = list(nx.get_node_attributes(G,"pos").values())

    raster_coords = make_LCZ_coords(raster_data,transform)

    # Use cKDTree for efficient nearest-neighbor search
    tree = KDTree(node_positions)
    distances, indices = tree.query(raster_coords)

    assigned_nodes = [list(G.nodes)[i] for i in indices]

    #  Create a new array with node assignments
    assigned_node_grid = np.reshape(assigned_nodes, raster_data.shape)

    # Filter by distance
    flat_data = raster_data.flatten()
    d_map = distances.copy()
    for i in range(len(distances)):
        dist = distances[i]
        node = assigned_nodes[i]
        lcz = flat_data[i]
        if dist <= max_dist: # If our node is too far away, we won't attribute demand to that location
            if "LCZ" not in G.nodes[node]:
                G.nodes[node]['LCZ'] = {}
            if lcz not in G.nodes[node]['LCZ']:
                G.nodes[node]['LCZ'][lcz] = 0
            G.nodes[node]['LCZ'][lcz] += 1 # Apply the LCZ info to the node
        else:
            d_map[i] = -1

    d_map = np.reshape(d_map, raster_data.shape)

    if show:
        raster_show(d_map,transform=transform,aspect=1,ax=ax) 
        ox.plot_graph(G,node_size=2,ax=ax)


# Given LCZ dict numbers make demand
def lcz_dict_to_weight(lcz_dict,avg_volume,total_volume):
    percent_of_the_pie = 0
    for lcz,num in lcz_dict.items():
        vol = avg_volume.get(lcz,0)
        percent_of_the_pie += num*vol/total_volume
    return percent_of_the_pie

def find_lcz_volume(lcz_dict,avg_volume):
    volume = 0
    for lcz,num in lcz_dict.items():
        
        vol = avg_volume.get(lcz,0)
        volume += num*vol
    return volume

def find_total_volume(G,avg_volume):
    total_volume = 0
    for lcz_dict in nx.get_node_attributes(G,"LCZ").values():
        total_volume+=find_lcz_volume(lcz_dict,avg_volume)
    return total_volume

def lcz_to_demand(G,avg_volume,total_demand):
    total_volume = find_total_volume(G,avg_volume)
    for node,lcz_dict in nx.get_node_attributes(G,"LCZ").items():
        G.nodes[node]['demand'] = lcz_dict_to_weight(lcz_dict,avg_volume,total_volume) * total_demand


def get_avg_volume():
    avg_volume = {1:12.5,2:9.625,3:1.95,4:7.5,5:5.25,6:4.875,7:1.2,8:0.975,9:1.625,10:25} # m^3/m Building volume to total volume From Rehm et al. 202
    for key,val in avg_volume.items():
        avg_volume[key] *= 100*100 # 100 m^2 blocks
    return avg_volume


def assign_LCZ_to_G(G,filename,total_demand):

    fig,ax = plt.subplots(figsize=(10,6))
    raster_data, transform, crs = load_lcz_raster(filename)
    find_closest_node_to_LCZ(G,raster_data,transform,ax=ax,max_dist=0.002) # max_dist has a large effect on output
    avg_volume = get_avg_volume() # From Rehm et al. Total volume of building per LCZ
    lcz_to_demand(G,avg_volume=avg_volume,total_demand=total_demand) # Total demand needs to be found
    
    