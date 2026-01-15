import wntr as wn
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import networkx as nx
import math


# Demand
def generate_demand_contour(G,ylim=None,show=True):
    network_graph = G
    # Extract node coordinates and demands
    node_positions = {node: (data['x'], data['y']) for node, data in network_graph.nodes(data=True) if 'x' in data and 'y' in data}
    node_demands = {node: network_graph.nodes[node].get('demand', 0) for node in network_graph.nodes()}

    # Prepare data for contour plot
    nodes = list(node_positions.keys())
    x = np.array([node_positions[node][0] for node in nodes])
    y = np.array([node_positions[node][1] for node in nodes])
    demands = np.array([node_demands[node] for node in nodes])

    # Define grid for contour plot
    grid_x, grid_y = np.meshgrid(np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100))

    # Interpolate demands over the grid
    grid_demands = griddata((x, y), demands, (grid_x, grid_y), method='linear')
    # print(demands)
    # Plot contour
    plt.figure(figsize=(8, 6))
    if not ylim:
        ylim = (min(demands),max(demands))
    plt.contourf(grid_x, grid_y, grid_demands, cmap='viridis',vmin=ylim[0],vmax=ylim[1])
    
    plt.colorbar(label='Demand')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Contour Plot of Node Demands')
    if show:
        nx.draw(network_graph, pos=node_positions, with_labels=False, node_color='red',
             node_size=3, edge_color='black', alpha=0.5,arrows=False)
        plt.show()


    return grid_x, grid_y, grid_demands,ylim


def linear_scale_coordinates(source_network, target_network):
    # Get minimum and maximum x and y values from both networks
    source_min_x = min(data['x'] for _, data in source_network.nodes(data=True) if 'x' in data)
    source_max_x = max(data['x'] for _, data in source_network.nodes(data=True) if 'x' in data)
    source_min_y = min(data['y'] for _, data in source_network.nodes(data=True) if 'y' in data)
    source_max_y = max(data['y'] for _, data in source_network.nodes(data=True) if 'y' in data)

    target_min_x = min(data['x'] for _, data in target_network.nodes(data=True) if 'x' in data)
    target_max_x = max(data['x'] for _, data in target_network.nodes(data=True) if 'x' in data)
    target_min_y = min(data['y'] for _, data in target_network.nodes(data=True) if 'y' in data)
    target_max_y = max(data['y'] for _, data in target_network.nodes(data=True) if 'y' in data)

    # Calculate scaling factors
    scale_x = (target_max_x - target_min_x) / (source_max_x - source_min_x)
    scale_y = (target_max_y - target_min_y) / (source_max_y - source_min_y)

    # Apply linear transformation to make source look like target network coordinates
    for node, data in source_network.nodes(data=True):
        if 'x' in data:
            data['x'] = target_min_x + scale_x * (data['x'] - source_min_x)
        if 'y' in data:
            data['y'] = target_min_y + scale_y * (data['y'] - source_min_y)

    return scale_x, scale_y



def apply_demand_contour_to_network(contour_x, contour_y, contour_demands, target_network,if_nan):
    # Interpolate contour demands onto new network nodes
    for node in target_network.nodes():
        x, y = target_network.nodes[node].get('x', 0), target_network.nodes[node].get('y', 0)
        if x != 0 and y != 0:  # Check if node has coordinates
            # Find the closest point on the contour grid
            closest_index = np.argmin((contour_x - x)**2 + (contour_y - y)**2)
            contour_demand = contour_demands.flatten()[closest_index]

            # Assign the closest contour demand to the node
            target_network.nodes[node]['demand'] = float(change_nan(contour_demand,if_nan))
        if target_network.nodes[node]['demand'] <= if_nan:
            target_network.nodes[node]['demand'] = if_nan

def change_nan(val, if_nan):
    if math.isnan(val):
        return if_nan
    else:
        return val

def attach_demands(G0,G1,if_nan):

    scale_x,scale_y = linear_scale_coordinates(G0, G1)

    contour_x, contour_y, contour_demands, ylim = generate_demand_contour(G0)



    # Assuming 'target_network' is the new network where you want to apply the demand contour
    apply_demand_contour_to_network(contour_x, contour_y, contour_demands, G1,if_nan)

    contour_x, contour_y, contour_demands, ylim = generate_demand_contour(G1)

# Features
def find_nearest_node(G, x, y):
    # Calculate distances from (x, y) to all nodes in the network
    distances = {node: np.sqrt((G.nodes[node]['x'] - x)**2 + (G.nodes[node]['y'] - y)**2)
                 for node in G.nodes() if 'x' in G.nodes[node] and 'y' in G.nodes[node]}

    # Find the node with the minimum distance
    nearest_node = min(distances, key=distances.get)
    min_distance = distances[nearest_node]

    return nearest_node

def find_distance(point1,point2):
    return ((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)**0.5

def copy_special_node_locations(G0,G1):
    """So reservoirs and tanks aren't overwritten"""
    scale_x,scale_y = linear_scale_coordinates(G0, G1)
    print(scale_x,scale_y)
    new_G1 = G1.copy()

    for res in nx.get_node_attributes(G0,"special"):
        if G0.nodes[res]['special'] == 'Reservoir':
            x = G0.nodes[res]['x']  
            y = G0.nodes[res]['y'] 
            nn = find_nearest_node(G1,x,y)


            attr = G0.nodes[res]
            if 'base_head' in attr:
                del attr['base_head']
            new_G1.add_node(res,elevation=G1.nodes[nn]['elevation'],**attr)


            new_G1.nodes[res]['x'] = x
            new_G1.nodes[res]['y'] = y
            new_G1.nodes[res]['demand'] = 0
            l = find_distance((x,y),(G1.nodes[nn]['x'],G1.nodes[nn]['y']))
            new_G1.add_edge(nn,res,length=l,pumphead=1,diameter=0.6,special='pump',pump=True)
        elif G0.nodes[res]['special'] == 'Tank':
            x = G0.nodes[res]['x']
            y = G0.nodes[res]['y']
            nn = find_nearest_node(G1,x,y)
            attr = G0.nodes[res]
            
            del attr['min_level']
            del attr['max_level']
            del attr['init_level']

            new_G1.add_node(res,**attr)
            new_G1.nodes[res]['x'] = x
            new_G1.nodes[res]['y'] = y
            new_G1.nodes[res]['demand'] = 0
            l = find_distance((x,y),(G1.nodes[nn]['x'],G1.nodes[nn]['y']))
            new_G1.add_edge(nn,res,length=l,diameter=0.6)
        else:
            print(G0.nodes[res])
    return new_G1
    
    