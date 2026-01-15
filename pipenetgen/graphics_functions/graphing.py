# graphics_functions/graphing.py
import osmnx as ox
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt



def plot_network(G,edge_width,edge_color,ax=None):
    """Plot the model network, # TODO: Move to graphics functions

    Args:
        G (osmnx graph): the graph to plot
        edge_width (str, optional): the parameter to graph on the edge. Defaults to None.
        edge_color (str, optional): the parameter to graph on the edge. Defaults to None.
    """
    nc = ["r" if G.nodes[node]['demand'] <= 0  else "k" for node in G.nodes] # Node colors are demand based
    ns = [abs(G.nodes[node]['demand']*4) if G.nodes[node]['demand'] >= 0 else 10 for node in G.nodes] 

    if edge_width:
        es = np.array(list(nx.get_edge_attributes(G,edge_width).values()))*10
    else:
        es = 1.5

    if edge_color:
        ecs = ox.plot.get_edge_colors_by_attr(G,edge_color,cmap='viridis_r')
    else:
        ecs = 'gray'
    if ax:
        ox.plot_graph(G,node_color=nc,node_size=ns,bgcolor=(0,0,0,0),edge_linewidth=es,edge_color=ecs,show=False,ax=ax)
        ax.set_aspect("equal")
        plt.show()
    else:
        ox.plot_graph(G,node_color=nc,node_size=ns,bgcolor=(0,0,0,0),edge_linewidth=es,edge_color=ecs,show=True,ax=ax)

import osmnx as ox
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.cm as mplcmap
from matplotlib.colors import Normalize
import numpy as np
from shapely.geometry import LineString
import matplotlib.patches as mpatches

def get_color_by_demand(demand, min_demand, max_demand):
    # If demand is negative, return gray
    if demand < 0:
        return 'gray'
    
    # Normalize demand to the range [0, 1]
    norm_demand = (demand - min_demand) / (max_demand - min_demand)
    
    # Get color from plasma colormap
    cmap = plt.cm.plasma
    color = cmap(norm_demand)
    
    # Convert RGBA to hex
    return color


def pretty_plot(G,ax=None,proj=True,mult=20,legend=True,node_mult=100,node_color=None,legend_shift=1.3):
    G = G.copy()
    nx.convert_node_labels_to_integers(G,0)

    if ax:
        show=False
    else:
        fig,ax = plt.subplots(figsize=(12,12))
        show=True
    if proj:
        G = ox.project_graph(G)  # Converts the graph to a projected CRS (UTM)
        

    # Define edge diameter bins and corresponding colors
    diameter_bins = [0.1016, 0.2032, 0.3048, 0.4064, 0.6]  # Adjust based on your data
    bin_indices = range(len(diameter_bins) - 1)

    # Normalize bins for colormap
    norm = Normalize(vmin=min(diameter_bins), vmax=max(diameter_bins))
    colormap = mplcmap.tab10

    # Prepare node attributes
    node_colors = []
    node_sizes = []
    node_alpha = []
    max_demand = max(list(nx.get_node_attributes(G,'demand').values()))
    for node, data in G.nodes(data=True):
        if data.get('special') == 'Reservoir':
            node_colors.append('blue')
            node_sizes.append(100)  # Fixed size for reservoirs
            node_alpha.append(0)
        elif data.get('special') == 'Tank':
            node_colors.append('red')
            node_sizes.append(100)  # Fixed size for tanks
            node_alpha.append(0)
        else:
            if node_color:
                node_colors.append(get_color_by_demand(data['demand'],0,max_demand))
            else:
                node_colors.append('black')
            node_sizes.append(data['demand'] * node_mult)  # Scale size by demand
            node_alpha.append(1)
    

    # Prepare edge attributes
    edge_colors = []
    edge_widths = []
    pump_edges = []  # Store pump edges for separate plotting

    for u, v, data in G.edges(data=True):
        if data.get('special') == 'pump':
            edge_colors.append((1,0,0,0.5))  # Pumps are red
            edge_widths.append(4)  # Thicker for pumps
            pump_edges.append((u, v))  # Collect pump edges

        else:
            diameter = data['diameter']
            # Find the appropriate color for the diameter
            bin_index = next((i for i, b in enumerate(diameter_bins) if diameter <= b), len(diameter_bins) - 1)
            color = colormap(norm(diameter_bins[bin_index]))  # Get the color from colormap
            
            edge_colors.append(color)
            edge_widths.append(diameter*mult)  # Normal edge width
   
    # Plot using OSMNX
    

    ox.plot_graph(
        G,
        node_color=node_colors,
        node_size=node_sizes,
        edge_color=edge_colors,
        edge_linewidth=edge_widths,
        node_alpha = node_alpha,
        bgcolor=(0,0,0,0),
        show=False,
        close=False,
        ax=ax
    )
    if not proj:
        ax.set_aspect(1)
    # Add legends
    # Node legend
    node_legend_handles = [
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='blue', markersize=10, label='Reservoir'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='red', markersize=10, label='Tank'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Demand Node'),
    ]
    # ax.legend(handles=node_legend_handles, loc='upper left', fontsize=10, title='Node Types')

    # Edge 
    from matplotlib.lines import Line2D

    edge_legend_handles = []
    # pump_handle_dashed = plt.Line2D([0], [0], color='red', lw=4, label='Pump',ls=':')
    pump_handle_solid = plt.Line2D([0], [0], color='red', lw=4, label='Pump',alpha=0.5)

    diameter_names = [4,8,12,16,24]

    pipe_legend_handles = [
        Line2D([0], [0], color=colormap(norm(diameter)), lw=diameter_bins[i]*mult, label=f"{diameter_names[i]} in")
        for i, diameter in enumerate(diameter_bins[:])
    ]
    edge_legend_handles.append(pump_handle_solid)
    edge_legend_handles.extend(pipe_legend_handles)


    # Add scale bar
    # Get the figure and axis for plotting

    from matplotlib_scalebar.scalebar import ScaleBar
    if proj:
    # Determine the bounds of the graph in projected coordinates
        x_min, x_max, y_min, y_max = ax.axis()

        # Calculate the distance between the min and max coordinates in meters (length of plot in x direction)
        scale_length = x_max - x_min
        print(scale_length)
        scalebar = ScaleBar(1, "m", location="upper left",length_fraction=0.2)


    # Symbology
    # Edges:

    # Overlay pump edges with dotted lines
    for u, v in pump_edges:
        x = [G.nodes[u]['x'], G.nodes[v]['x']]
        y = [G.nodes[u]['y'], G.nodes[v]['y']]
        ax.plot(x, y, color='red', linestyle=':', linewidth=4,zorder=0, label='Pump' if 'Pump' not in [l.get_label() for l in ax.get_legend_handles_labels()[0]] else None)

    # Nodes
    node_positions = {node: (data['x'], data['y']) for node, data in G.nodes(data=True)}

    reservoir_nodes = [node for node, data in G.nodes(data=True) if data.get('special') == 'Reservoir']
    tank_nodes = [node for node, data in G.nodes(data=True) if data.get('special') == 'Tank']

    reservoir_positions = [node_positions[node] for node in reservoir_nodes]
    ax.scatter(
        [pos[0] for pos in reservoir_positions],
        [pos[1] for pos in reservoir_positions],
        c='blue', s=100, marker='s', label='Reservoir',edgecolors='k'
    )

    # Tanks
    tank_positions = [node_positions[node] for node in tank_nodes]
    ax.scatter(
        [pos[0] for pos in tank_positions],
        [pos[1] for pos in tank_positions],
        c='red', s=100, marker='^', label='Tank',edgecolors='k'
    )

    # Node legend
    import matplotlib.lines as mlines

    reservoir_legend = mlines.Line2D([],[],marker='s',color='blue', label='Reservoir', markeredgecolor='black', linewidth=0,markersize=10)
    tank_legend = mlines.Line2D([],[],markersize=10,color='red', marker='^',label='Tank', markeredgecolor='black', linewidth=0)
    node_legend = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=8, label='Node')

    handles =[reservoir_legend, tank_legend, node_legend,pump_handle_solid, mlines.Line2D([],[],linewidth=0)]
    handles.extend(pipe_legend_handles)
    print(handles)
    # Add the legend to the plot 
    if legend:
        ax.legend(
            handles=handles,
            loc='upper right',
            fontsize=10,
            ncol = 2, 
            frameon = True, 
            bbox_to_anchor = (legend_shift,1),
            prop={'size':20}
        )


    # Show the plot
    if show:
        plt.show()
    
def plot_error_cost(errors,costs):
    # Example data
    costs = [c/1000000 for c in costs]  # Replace with actual cost values
    
    fig, ax1 = plt.subplots(figsize=(6,3))
    costs = costs[1:]
    errors = errors[1:] 
    iterations = list(range(1, len(costs) + 1))  # Assuming costs and errors are the same length


    # Plot costs on the primary y-axis
    ax1.plot(iterations, costs, 'C0-o', label='Cost')  # Black line with circles
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Cost (M$)', color='k')
    ax1.tick_params(axis='y', colors='k')

    # Create secondary y-axis for errors
    ax2 = ax1.twinx()
    ax2.plot(iterations, errors, 'C1--s', label='Error')  # Black dashed line with squares
    ax2.set_ylabel('Error (m)', color='k')
    ax2.tick_params(axis='y', colors='k')

    # Create a combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

    fig.tight_layout()
    plt.xticks(range(1,len(iterations)+1))

    plt.show()