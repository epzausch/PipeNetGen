import pipenetgen as png
import osmnx as ox
import networkx as nx
import pandas as pd
import shapely
import time

                # max and min pressure
parameters = {"MinP": 20, # m
              "MaxP":56, # m
              # base demand
              "demand":0.02, # Lps average demand
              "elev_mult": 1, # elevation multiplier
              # Diam and costs
              "diameters":[0.1016, 0.2032, 0.3048, 0.4064, 0.6], # m
              "cost":[200,219,300,400, 634], # $/m
                # 'cost':[200,220,240,260,280],
              # default velocity
              "velocity":1, # m/s
              # EPANET ITerations
              'iterations':5,
              # FIRe FLow
              "hydr_demand":30,
              # Max and min flow
              "QMAX":10, # LPS, must be more than max flow reservoir provides
              "QMIN":0.1, # LPS
              # TOPO
              "link_density":0.0008,  #0.0008
              "node_degree":2, # 2
              "meshedness":0.1, #0.04
              # BIG M
              "BigM_P":200,
              "BigM_Q":11,
              # Pipe roughness
              "C":120, # Roughness Coeff
              # Tower costs
              "max_tower":30,
              "min_tower":10,
              "TowerCostPerm":100000,
              "TowerCost":1000000,
              # MILP params
              "Gap%":.1,
              "TimeLimit":600,
              "MIPFocus":1,
              "verbose":True,
              # Special features
              "min_velo":2,
              'all_edges':True,
              "max_head":85,
              'seed':0
              }


cluster_params = {"n_clusters":6,
                  "QMIN":0.002,
                  "QMAX":18,
                  "BigM_Q":19,
                  "Gap%":0.24,
                  "min_velo":1,
                  "all_edges": False,
                  "verbose": True,
                  "max_head":100}

fire_params = {'verbose': True,
               'QMAX': 30,
               'BigM_Q': 31,
               'Gap%': 0.01,
               'TimeLimit': 60*60,
               'min_velo': 1}

api_key = 'INSERT HERE IF USING IUWM DEMANDS'

inputs = {"polygon":{"polygon_filepath":"Input/LakewoodCoords.csv",'City_name':"Lakewood, California"},
          "tower":{"tower_filepath":"Input/lakewood_water_towers.csv"},
          "source":{'source_filepath':"Input/Lakewood_sources.csv",'num_source':8},
          "demand":{"demand_filepath":'Input/Lakewood_lcz.tif','type':"LCZ",'total_demand':18}}

save_to = {"hnt":"Output/Lakewood_folder/png_Lakewood_Hint_{0}.hnt",
           "fire":"Output/Lakewood_folder/png_testing_hnt_lakewood_fire.hnt",
           'final':"Output/Lakewood_folder/png_Final_Lakewood"}

iteration_params = {'QMIN':0}


def pipe_network_generator(inputs,save_to,parameters=parameters,cluster_params=cluster_params,iteration_params=iteration_params,fire_params=fire_params):
    """Generates a pipe network from inputs and saves as an inp

    Args:
        inputs (_type_): _description_
        save_to (_type_): _description_
        parameters (_type_, optional): _description_. Defaults to parameters.
    """
    parameters = parameters.copy()
    cluster_params = cluster_params.copy()
    fire_params= fire_params.copy()

    start_time = time.time()
    # Create the network and prepare it
    G = network_prep(inputs,parameters)
    png.plot_network(G,'length','length')
    # Cluster the network
    G = cluster_the_network(G,parameters,cluster_params)

    if inputs['source'].get('real') == True:
        G = update_attrs_for_real(G,parameters)

    G = nx.convert_node_labels_to_integers(G)


    
    # MILP - EPANET Iteration
    G = MILP_EPANET(G,save_to.get("hnt"),iteration_params,parameters)
    # Fire Flow
    G = fire_flow(G,save_to.get("hnt"),save_to.get("fire"),parameters,fire_params)

    end_time = time.time()
    elapsed_time = end_time - start_time
    G.graph['time'] = elapsed_time

    # Save
    save_network(G,save_to.get('final'))

    return G

def update_attrs_for_real(G,parameters):
    for edge in nx.get_edge_attributes(G,'special'):
        # print(G2.edges[edge])
        u,v,k = edge
        G.edges[edge]['pumphead'] = parameters['max_head'] - G.nodes[u]['elevation']

    for node in nx.get_node_attributes(G,'special'):
        if "vol_curve_name" in G.nodes[node]:
            del G.nodes[node]['vol_curve_name']
            print("del")
        G.nodes[node]['head'] = G.nodes[node]['elevation']
        G.nodes[node]['pressure'] = 0
    return G


def save_network(G,final_save):
    png.G_to_inp(G,final_save+".inp")
    png.G_to_osm(G,final_save+".osm")

import networkx as nx

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


def fire_flow(G,hint_read_write,fire_hint_path,parameters,fire_params):
    
    hydr_locations = list(png.find_hydrant_locations(G,fire_params['fire_gap']))
    try:
        hydr_locations = sort_by_proximity_to_special(G.to_undirected(),hydr_locations)
    except:
        print("hydrant nodes not sorted")
    # Highlight hydrant locations
    nc = ['r' if node in hydr_locations else "w" for node in G.nodes]
    ox.plot_graph(G, node_color=nc, node_size=50, node_zorder=2, show=True)
    for node in hydr_locations:
        G.nodes[node]['hydrant'] = True
    for node in G.nodes:    
        G.nodes[node]['base_demand'] = G.nodes[node]['demand']
    # G.nodes[hydr_locations[0]]['demand'] += hydr_demand
    png.set_reservoir_pressure_to_pump(G)
    epa_G_ff = png.run_epanet_algorithm(G,parameters)

    parameters.update(fire_params)

    errors_ff = [png.calculate_pressure_dif(epa_G_ff)]
    costs_ff = [png.calc_cost(epa_G_ff,parameters)]

    ct=  0
    for i in range(len(hydr_locations)-1):
        if i == 0:
            hnt_path = hint_read_write.format(parameters.get('iterations'))
        else:
            hnt_path = fire_hint_path
        try:
            epa_G_ff,errors_ff,costs_ff = png.run_a_fire_flow(epa_G_ff,parameters,hydr_locations[i],hnt_path,fire_hint_path,errors_ff,costs_ff,False)
            png.plot_network(epa_G_ff,'diameter','flowrate')
            epa_G_ff.nodes[hydr_locations[i]]['hydrant_ran'] = True
        except:
            epa_G_ff.nodes[hydr_locations[i]]['hydrant_ran'] = False
            print("MODEL FAILED\n\n\n\n")
            print(i,hydr_locations[i])
            epa_G_ff.nodes[hydr_locations[i]]['demand'] = epa_G_ff.nodes[hydr_locations[i]]['base_demand']
            ct += 1
            print("\n\n\n\nMODEL FAILED")

    for node,dem in nx.get_node_attributes(epa_G_ff,'demand').items():
        if dem >= 30:
            epa_G_ff.nodes[node]['demand'] = epa_G_ff.nodes[node]['base_demand']


    png.plot_error_cost(errors_ff,costs_ff)
    return epa_G_ff

def MILP_EPANET(G,hint_read_write,iteration_params,parameters):
    # Run EPANET once
    parameters.update(iteration_params)

    if G.graph['cluster'] == True:
        epa_G = png.run_epanet_algorithm(G,parameters)
        errors = [0]
        costs = [0]
    else:
        epa_G = G
        errors = [0]
        costs = [png.calc_cost(epa_G,parameters)]

    
    max_iterations = parameters.get("iterations", 5) # Never needs to go to 5
    min_change_ratio = 0.10  # 10%

    for i in range(1, max_iterations + 1):

        # Determine whether we should write the HNT file this iteration
        if G.graph['cluster'] == False:
            write_hnt = (i == 2)
        else:
            write_hnt = (i == 1)

        # Run EPANET step
        epa_G, errors, costs = png.run_a_run(
            epa_G,
            parameters,
            hint_read_write.format(i - 1),
            hint_read_write.format(i),
            errors,
            costs,
            write_hnt
        )

        png.plot_network(epa_G, 'diameter', 'flowrate')

        # ---- STOPPING CONDITION ----
        # Need at least two errors to compare
        if len(errors) >= 2:
            prev_err = errors[-2]
            curr_err = errors[-1]

            # Avoid division-by-zero
            if prev_err != 0:
                change_ratio = abs(curr_err - prev_err) / abs(prev_err)

                # Stop if error changed by less than 10%
                if change_ratio < min_change_ratio:
                    print(f"Stopping early at iteration {i}: error change {change_ratio:.4f}")
                    break

        # Also stop if we reached max iterations
        if i == max_iterations:
            print(f"Reached max iterations ({max_iterations}).")
            break

    png.plot_error_cost(errors,costs)
    return epa_G





def cluster_the_network(G,parameters,cluster_params):
    if cluster_params['n_clusters'] > 1:
        base = png.run_cluster_algorithm(G,cluster_params['n_clusters'],parameters,show=cluster_params['verbose'],cluster_params=cluster_params,target_topo=parameters)
        png.plot_network(base,"diameter",'flowrate')
        base.graph['cluster'] = True
        return base
    else:
        G.graph['cluster'] = False
        for edge in list(G.edges):
            u,v,k = edge
            if k != 0:
                G.add_edge(u,v,0,**G.edges[edge])
                G.remove_edge(u,v,k)
        return G

def network_prep(inputs,parameters):
    """Creates and prepares a network for optmization based on user defined inputs

    Args:
        inputs (_type_): _description_
        parameters (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Load Network
    G = load_streets(inputs['polygon'].get('polygon_filepath'),inputs['polygon'].get("City_name")) # TODO: Raise Custom error if nothing found
    

    # If we don't have source locations we need to do this info twice.
    # png.assign_meteo_elevation(G,grp_size=10)
    png.assign_open_elevation(G)

    # Load Sources
    updated = False
    if inputs['source'].get("source_filepath",False) == False and inputs['source'].get('real',False) == False:
        print("updating early")
        update_topo(G)
        updated = True

    G = load_sources(G,inputs['source'].get("source_filepath"),inputs['source'].get('num_source'),inputs['source'])

    # Load Towers
    if inputs.get('tower'):
        png.add_tanks_from_csv(G,inputs['tower'].get('tower_filepath'))
    ox.plot_graph(G)
    if updated == False:
        # Fix topology
        G = update_topo(G)

    # Add elevation
    # png.assign_meteo_elevation(G,grp_size=10)
    png.assign_open_elevation(G)
    # Add Demand
    G = give_demands(G,inputs['demand'].get("demand_filepath"),inputs['demand'].get("type"),inputs['demand'].get('total_demand'),parameters['demand'],inputs['demand'])

    

    G.graph['max_head'] = parameters['max_head']
    # print("MAX HEAD",G.graph['max_head'])

    return G.to_undirected()


def give_demands(G,demand_filepath,dem_type,total_demand,demand,demand_params):
    png.add_demands(G,demand)
    if demand_params.get('real') == True:
        png.attach_demands(png.load_real(),G,if_nan=0.2)
    else:
        if demand_filepath:
            if dem_type == "LCZ":
                png.assign_LCZ_to_G(G,demand_filepath,total_demand)
            elif dem_type == "IUWM":
                api_key = parameters.get("api_key")
                png.attach_IUWM_demands(G,demand_filepath,api_key)
        else:
            if total_demand:
                set_demand = total_demand/len(G.nodes)
            else:
                set_demand = demand
            for node in G.nodes:
                G.nodes[node]['demand'] = set_demand
        for node in G.nodes:
            if G.nodes[node]['demand'] <= demand and G.nodes[node]['demand'] >= -0.0001:
                G.nodes[node]['demand'] = demand
    
    print("DEMANDS")
    print(sum([G.nodes[node]['demand'] for node in G.nodes  if G.nodes[node]['demand'] > 0 ]))
    png.redistribute_demands(G)
    if total_demand:
          current_demand = sum([G.nodes[node]['demand'] for node in G.nodes  if G.nodes[node]['demand'] > 0 ])
          multiplier = total_demand/current_demand
          for node in G.nodes:
              G.nodes[node]['demand'] *= multiplier
    print(sum([G.nodes[node]['demand'] for node in G.nodes  if G.nodes[node]['demand'] > 0 ]))
    print(sum([G.nodes[node]['demand'] for node in G.nodes]))

    return G


def update_topo(G):
    G = nx.convert_node_labels_to_integers(G)
    G = png.fix_topo(G,show=True,consolidate_dist=15)
    return G

def load_sources(G,source_filepath,n_sources,source_params):
    if source_filepath:
        png.add_sources_from_csv(G,source_filepath)
    else:
        if source_params.get('real') == True:
            G = png.copy_special_node_locations(png.load_real(),G)
        else:
            png.add_sources(G,n=n_sources)
    return G

def load_streets(polygon_filepath,City_name):
    # Load polygon
    if polygon_filepath:
        coords = pd.read_csv(polygon_filepath)
        poly = shapely.geometry.Polygon(coords[['x','y']])
        # Create Graph
        G = ox.graph_from_polygon(poly,network_type="drive")
    else:
        G = ox.graph_from_place(City_name,network_type='drive')

    return G
