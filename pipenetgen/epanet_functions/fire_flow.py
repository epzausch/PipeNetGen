# Functions for running fire flow simulation, no longer the recomended method (see full_demo).
import wntr
import osmnx as ox
import networkx as nx


# Function to find the optimal hydrant locations
def find_hydrant_locations(G, max_distance):
    # Convert the graph to undirected and get the adjacency list
    H = G.to_undirected()
    adj = {n: set(H.neighbors(n)) for n in H.nodes()}
    
    # Set of nodes to cover
    nodes_to_cover = set(H.nodes())-set(nx.get_node_attributes(G,'special').keys())
    
    # Set to store chosen hydrant locations
    hydrants = set()
    
    while nodes_to_cover:
        # Find the node with the maximum neighbors within max_distance
        best_node = None
        best_coverage = 0
        for node in nodes_to_cover:
            # Use NetworkX to get nodes within max_distance edges
            coverage = set(nx.single_source_shortest_path_length(H, node, cutoff=max_distance).keys())
            if len(coverage & nodes_to_cover) > best_coverage:
                best_node = node
                best_coverage = len(coverage & nodes_to_cover)
        
        # Add the best node to hydrants
        hydrants.add(best_node)
        
        # Remove covered nodes from nodes_to_cover
        coverage = set(nx.single_source_shortest_path_length(H, best_node, cutoff=max_distance).keys())
        nodes_to_cover -= coverage
    
    return hydrants



def sort_by_proximity_to_special(graph, hydr_locations):
    """
    Sort nodes based on proximity to the nearest node with a 'special' attribute.

    Parameters:
        graph (networkx.Graph): The graph containing nodes.
        hydr_locations (list): List of nodes to sort.

    Returns:
        list: Sorted list of nodes by proximity to the nearest 'special' node.
    """
    # Step 1: Find all nodes with the 'special' attribute
    special_nodes = [n for n, attrs in graph.nodes(data=True) if attrs.get('special')]

    if not special_nodes:
        raise ValueError("No nodes with the 'special' attribute found in the graph.")

    # Step 2: Compute shortest distances to the closest 'special' node
    def distance_to_special(node):
        return min(nx.shortest_path_length(graph, source=node, target=sp) for sp in special_nodes)

    # Step 3: Sort the hydr_locations based on distances to special nodes
    sorted_hydr_locations = sorted(hydr_locations, key=distance_to_special)
    
    return sorted_hydr_locations
