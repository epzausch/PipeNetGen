# Functions for running wntr simulations
import wntr
import osmnx as ox
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
# from ..data_functions import G_to_wn
import math
from ..preperation_functions import *
from ..linear_programs import *
import shapely as sp
from gurobipy import GRB

def delete_pump_head(G):
    for edge in nx.get_edge_attributes(G,"pumphead"):
        if 'pumphead' in G.edges[edge]:
            del G.edges[edge]['pumphead']
        if 'pump' in G.edges[edge]:
            del G.edges[edge]['pump']
        if 'special' in G.edges[edge]:
            del G.edges[edge]['special']

def run_a_fire_flow(epa_G,parameters,hydrant_location,hnt_path,save_to,errors,costs,write=False):
    epa_G.nodes[hydrant_location]['demand'] += parameters['hydrant_demand']
    epa_G = run_epanet_algorithm(epa_G,parameters)
    delete_pump_head(epa_G)
    # RUN LP
    # png.node_elev_mult(epa_G,parameters['elev_mult'])
    lp = milp(epa_G,parameters,network_type="fire")

    seed = parameters.get('seed',0)
    lp.m.params.seed = seed

    lp.load_model()

    # hnt_path = f'Output/determinism/{seed_n}_{iter_n}_last_solution.hnt'
    # png.write_hnt_file(lp,hnt_path)
    if write == True:
        write_hnt_file(lp,hnt_path)

    lp.m.read(hnt_path)
    

    lp.run_model()
    lp.m.update()
   
    costs.append(lp.m.getObjective().getValue())
    # lp.m.setParam("LogFile", f"Output/determinism/{seed_n}_log.log")

    if lp.m.status == GRB.OPTIMAL:
        for v in lp.m.getVars():
            v.VarHintVal = v.X
            v.VarHintPri = 1
        
        
        lp.m.write(save_to)

    
    lp.G.nodes[hydrant_location]['demand'] = lp.G.nodes[hydrant_location]['base_demand']# parameters['hydrant_demand']
    

    epa_G2 = run_epanet_algorithm(lp.G,parameters)

    # SAVE
    errors.append(calculate_pressure_dif(epa_G2))
    
    return epa_G2,errors,costs

   

def run_a_run(epa_G,parameters,hnt_path,save_to,errors,costs,write=True):

    delete_pump_head(epa_G)
    # RUN LP
    # png.node_elev_mult(epa_G,parameters['elev_mult'])
    lp = milp(epa_G,parameters,network_type="simple")

    seed = parameters.get('seed',0)
    lp.m.params.seed = seed

    lp.load_model()
    

    # hnt_path = f'Output/determinism/{seed_n}_{iter_n}_last_solution.hnt'
    # png.write_hnt_file(lp,hnt_path)
    if write == True:
        write_hnt_file(lp,hnt_path)
    else:
        try:
            lp.m.read(f"{hnt_path}.sol")
            lp.m.read(f"{hnt_path}.mst")
        except:
            pass
    try: 
        lp.m.read(hnt_path)
    except:
        pass
    
    lp.run_model()
    lp.m.update()
    if write == True:
        lp.m.write(f"{save_to}.sol")
        lp.m.write(f"{save_to}.mst")
        lp.m.write(save_to)

    costs.append(lp.m.getObjective().getValue())
    # lp.m.setParam("LogFile", f"Output/determinism/{seed_n}_log.log")

    # if lp.m.status == GRB.OPTIMAL:
    for v in lp.m.getVars():
        v.VarHintVal = v.X
        v.VarHintPri = 1
    
    
    epa_G2 = run_epanet_algorithm(lp.G,parameters)
    #epa_G2 = lp.G
    # SAVE
    errors.append(calculate_pressure_dif(epa_G2))
    
    return epa_G2,errors,costs

def G_to_wn(G,geometry=True):
    """Convert a osmnx graph into a wntr network

    Args:
        G (osmnx graph): an osmnx graph
        geometry (bool, optional): Whether you want to preserve geometry. Defaults to True.

    Returns:
        wntr network: the same osmnx graph but in wntr, ready to be run in epanet
    """
    if G.graph.get('max_head'):
        max_head = G.graph['max_head']
        # print(max_head)
        # print("Max Head Set")
    elif len(nx.get_node_attributes(G,'head')) > 0:
        max_head = max(list(nx.get_node_attributes(G,'head').values()))
        # print("Max Head Calcualted")
    else: 
        # print("Missing Max Head")
        max_head = 85
    
    # Replace all nodes with ints

    wn = wntr.network.WaterNetworkModel()
    wn.add_pattern("base",[1]) 
    # Add nodes
    for node, data in G.nodes(data=True):
        spesh = G.nodes[node].get("special",0)
        if spesh in ['Reservoir','source']:
            # base_head = base_head=data.get('head', 50)
            # print(data.get('pressure', 0)+data.get('elevation',0))
            # print(data)
            wn.add_reservoir(str(node), base_head=data.get('pressure', 0)+data.get('elevation',0)) 
            joi = wn.get_node(str(node))
            
        elif spesh in ['Tank','tank','tower']: 
            # print(max_head,type(data.get('elevation')))
            max_head = float(max_head)

            init_level = data.get('init_level', max_head-float(data.get('elevation',0)))
            
            wn.add_tank(str(node), 
                        elevation=data.get('elevation', 0), 
                        init_level=init_level,
                        min_level=data.get('min_level', max_head-data.get('elevation',0)-10), # Tank is 10 m off the ground 
                        max_level=data.get('max_level', max_head-data.get('elevation',0)),
                        diameter=data.get('diameter', 10),
                        min_vol=data.get('min_vol', 0))
            joi = wn.get_node(str(node))
        else:
            wn.add_junction(str(node), base_demand=data.get('demand', 0)/1000,demand_pattern="base", demand_category="base")
            joi = wn.get_node(str(node))
        
        for k,v in data.items():
            
            if k not in ['base_head']:
                try:
                    setattr(joi,k, v)
                except:
                    pass
        # print(data['x'],data['lon'])
        joi.coordinates = [data["x"],data["y"]] 
       
       
    # Add edges
    for u, v, key, data in G.edges(keys=True, data=True):
        
        if 'special' in data:
            if data['special'] == 'pump' and not (G.nodes[u].get('special') and G.nodes[v].get('special')):
                pump_id = f'pump_{str(u)}_{str(v)}_{str(key)}'
                pump_name = pump_id + '_curve' 
                if len(pump_name) > 32:
                    pump_name = pump_name[:32] 
                # print("ADDING PUMP", pump_name)
                wn.add_curve(pump_name, 'HEAD', [(0, data.get('pumphead', 1)+0.0001), (10000, data.get('pumphead', 1))])
                wn.add_pump(pump_id, str(u), str(v), pump_parameter=pump_name, pump_type="HEAD")
                eoi = wn.get_link(f'pump_{u}_{v}_{key}')
            else:
                wn.add_pipe(f'pipe_{u}_{v}_{key}', str(u), str(v))
                eoi = wn.get_link(f'pipe_{u}_{v}_{key}')
        else:
            wn.add_pipe(f'pipe_{u}_{v}_{key}', str(u), str(v))
            eoi = wn.get_link(f'pipe_{u}_{v}_{key}')
        
        

        
        for key,val in data.items():
            # if k not in ["name","velocity"]:
            try:
                setattr(eoi,key, val)
            except:
                pass
        if geometry:
            vertices = list(data.get("geometry",calculate_geometry(G,u,v)).coords) 
            eoi.vertices = fix_vertices(vertices,(G.nodes[u]['x'],G.nodes[u]['y']),(G.nodes[v]['x'],G.nodes[v]['y']))
        
        
    wn.options.hydraulic.inpfile_units = "LPS"
    return wn

def fix_vertices(vertices, start, end):
    from math import sqrt

    def distance(p1, p2):
        return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    # Unpack start and end coordinates
    sx, sy = start
    ex, ey = end

    # Check distances
    dist_start = distance(vertices[0], (sx, sy))  # Distance from first vertex to start
    dist_end = distance(vertices[0], (ex, ey))  # Distance from last vertex to end
    
    # If the last vertex is "jumping" back to start or not close to end, reverse vertices
    if dist_end < dist_start:
        vertices.reverse()
    
    return vertices


def calculate_geometry(G,node1,node2):
    """ Calcualtes straight geometery"""
    p1 = sp.Point(G.nodes[node1]['x'],G.nodes[node1]['y'])
    p2 = sp.Point(G.nodes[node2]['x'],G.nodes[node2]['y'])
    return sp.LineString([p1,p2])

def round_demands(G,rounding=[6,6]):
    """Round demands from a graph without changing total demand
    Used by transfer flowrate because EPANET is inconsistent.

    Args:
        G (osmnx graph): The graph whose source demands you wnat to round
        rounding (list, optional): first is the actual source second is the error (removes float errors). Defaults to [4,8].
    """
    total = get_total_source_demands(G)
    for node in nx.get_node_attributes(G,"special"):
        G.nodes[node]['demand'] = round(G.nodes[node]['demand'],rounding[0])
    new_total = get_total_source_demands(G)
    dif = round(total-new_total,rounding[1])
    # print(dif)
    for node in nx.get_node_attributes(G,"special"):
        G.nodes[node]['demand'] += dif
        break


def transfer_flowrates(wn,G,rounding=[6,6]):
    """Transfer flowrates, pressure, and demand information from epanet simulations onto a osmnx graph
    Similar to wntr_onto_osmnx but actually works. 

    Args:
        wn (wntr network): a wntr network that you want to run
        G (osmnx graph): a osmnx graph to receive the wntr information

    Returns:
        osmnx graph: a new osmnx graph with hydraulic information
    """
    epa_G = G.copy()

    # Simulate hydraulics using the WNTRSimulator
    sim = wntr.sim.WNTRSimulator(wn)
    results_WNTR = sim.run_sim()
    
    # Save the information
    flowrates = results_WNTR.link['flowrate']
    pressures = results_WNTR.node['pressure']

    inflow = results_WNTR.node['demand']
    node_demand = results_WNTR.node['demand']
    # Transfer the flowrates 
    
    for edge in epa_G.edges:
        u,v,k = edge
        base = str(u)+"_"+str(v)+"_"+str(k)
        if "pump_" + base in flowrates:
            # pass
            epa_G.edges[edge]['flowrate'] = round(float(flowrates["pump_"+base].iloc[0]) * 1000,4)
        elif "pipe_"+base in flowrates:
            # print(epa_G.edges[edge])
            # print(flowrates['pipe_'+base])
            epa_G.edges[edge]['flowrate'] = round(float(flowrates["pipe_"+base].iloc[0]) * 1000,4) # LPS

        else:
            if epa_G.edges[edge].get('name') in flowrates:
                epa_G.edges[edge]['flowrate'] = round(float(flowrates[epa_G.edges[edge].get('name')].iloc[0]) * 1000,4)
            else:
                print("No Edge",base)

        if 'diameter' in epa_G.edges[edge]:
            epa_G.edges[edge]['velocity'] = calculate_velocity(epa_G.edges[edge]['diameter']*1000,epa_G.edges[edge]['flowrate']/1000)


        
    for node in epa_G.nodes:
        epa_G.nodes[node]['pressure'] = round(float(pressures[str(node)].iloc[0]),4)
        if 'special' in epa_G.nodes[node]:
            epa_G.nodes[node]['demand'] = inflow[str(node)].iloc[0] * 1000
        else:
            epa_G.nodes[node]['demand'] = node_demand[str(node)].iloc[0] * 1000
    
    round_demands(epa_G,rounding) # Round demands to 6 decimal places

    return epa_G

def  find_total_wntr_node_demand(noi):
    """Calculate the total demand within a wntr node. No longer used as part of the methodology

    Args:
        noi (a wntr node object): the node of interest (noi)

    Returns:
        float: the total demand 
    """
    total_demand = 0
    # print(noi)
    for pattern in noi.demand_timeseries_list:
        # print(pattern)
        total_demand += pattern.base_value
    return total_demand


def calculate_pressure_dif(epa_G,char="pressure_dif"):
    """Calcualte the absolute sum of a nodal attribute in the osmnx network epa_G 
    (usually for pressuredifference)

    Args:
        epa_G (osmnx graph): the osmnx graph with nodes which have attribute char
        char (str, optional): The charecter to absolute sum. Defaults to "pressure_dif".

    Returns:
        float: the total magnitude of char
    """
    t_error = 0
    for val in nx.get_node_attributes(epa_G,char).values():
        t_error += abs(val)
    return t_error

def attach_pressure_dif(milp_G,epa_G):
    """Calculate the pressure difference between two osmnx graphs and attach them to the second one

    Args:
        milp_G (osmnx graph): a first osmnx graph with nodal pressures
        epa_G (osmnx graph): the second omsnx graph with nodal pressures. Will receive the difference between the two
    """
    for node in milp_G.nodes:
        if 'special' not in milp_G.nodes:
            milp_pressure = milp_G.nodes[node]['pressure']
            epa_pressure = epa_G.nodes[node]['pressure']
            dif = milp_pressure-epa_pressure
        else:
            dif = 0
        
        epa_G.nodes[node]['pressure_dif'] = dif

def remove_peaking(net,peaking_factor=3):
    for node in net.nodes:
        net.nodes[node]['peak_demand'] = net.nodes[node]['demand']
        net.nodes[node]['demand'] = net.nodes[node]['demand']/net.graph.get('peaking_factor',peaking_factor)

    epa_G = run_epanet_algorithm(net,{'max_head':net.graph.get('max_head',105)})
    for node in epa_G.nodes:
        epa_G.nodes[node]['head'] = epa_G.nodes[node]['pressure'] + epa_G.nodes[node]['elevation']
    return epa_G

def run_epanet_algorithm(G,parameters):
    G = G.copy()
    if 'max_head' not in parameters:
        max_head = parameters['MaxP'] + max(G.nodes[node]['elevation'] for node in nx.get_node_attributes(G,"special")) #max(G.nodes[node]['head'] for node in nx.get_node_attributes(G,"special"))
        parameters['max_head'] = max_head
    else:
        max_head = parameters['max_head']

    max_head = float(max_head)
    # Adjust pressure for all source nodes to ensure they have the same head
    for node in nx.get_node_attributes(G,"special"):
        # Calculate the pressure needed to reach the max head
        G.nodes[node]['pressure'] = max_head - G.nodes[node]['elevation']
        # Set the new head using the updated pressure
        G.nodes[node]['head'] = G.nodes[node]['pressure'] + G.nodes[node]['elevation']

    correct_pressures(G,120)

    set_reservoir_pressure_to_pump(G)  


    wn = G_to_wn(G)
    epa_G = transfer_flowrates(wn,G)


    attach_pressure_dif(G,epa_G)
    # make_boxplot(epa_G)

    return epa_G

    # This does not change the output randomly 6062 error

def write_hnt_file(lp,filepath):
    model = lp.m
    add_edgeloss(lp.G)
    all_vars = model.getVars()
    hint_lines = []
    i = 0
    # Try starting with large pipes? 
    for var in all_vars:
        val = 0
        name = var.VarName
        if name[:4] == 'diam':
            s,u,v = name.split("_")
            s = 0.6

            if (int(u),int(v),0) in lp.G.edges:
                eoi = lp.G.edges[(int(u),int(v),0)]
                actual_diam = eoi['diameter']

                if s == actual_diam:
                    val = 1
                else:
                    val = 0

        
        
    hint_lines.insert(0,"# MIP variable hints")

    
    with open(filepath,'w') as file:
        for line in hint_lines:
            file.write(line)
            file.write("\n")



def add_edgeloss(epa_G):
    for edge in epa_G.edges:
        d = epa_G.edges[edge]['diameter'] * 1000
        C = epa_G.edges[edge].get('roughness',120)
        L = epa_G.edges[edge]['length']
        Q = epa_G.edges[edge]['flowrate']

        V = calculate_velocity(d,Q)
        hl = calculate_edgeloss(d,C,L,V)

        # if 'edgeloss' in epa_G.edges[edge]:
        #     print(hl,epa_G.edges[edge]['edgeloss'])
        epa_G.edges[edge]['edgeloss'] = hl

def calculate_edgeloss(diameter,c,length,velocity):
    """Takes a diameter in mm a roughness coeffecient and length in m to calculate headloss with H-W


    Args:
        diameter (float): the diameter of the pipe (mm)
        c (float): the roughness (H-W)
        length (float): the length of the pipe (m)
        velocity (float): the velocity in the pipe (m/s)

    Returns:
        float: the head loss of the edge
    """

    head_loss = 10.67 * (velocity*(diameter/2)**2*math.pi)**1.852 * c**(-1.852) * (diameter)**-4.8704 * length
    return head_loss

def calculate_velocity(diameter,flowrate):
    """Calculate velocity (m/s) from diameter (mm) and flowrate (m^3/s)

    Args:
        diameter (float): the diameter of interest (mm)
        flowrate (flow): the flowrate of interest (m^3/s)

    Returns:
        float: a velocity (m/s)
    """
    velo = flowrate / (math.pi * (diameter / 2000) ** 2)
    return abs(velo)

def correct_pressures(G,c=120):
    """Correct the pressures of the MILP output. Done roughly, not very accurate, room for marginal improvements.
    Why does the MILP need pressure correction? 

    Because pressure constraints between nodes only ensure that the second node is less than what it should be
    This means the MILP underestimates nodal pressure when the flowrate is correct.

    Args:
        G (osmnx graph): an osmnx graph run through the milp
        c (int, optional): the roughness (H-W). Defaults to 120.
    """
    dif = 10 # any non 0 number
    breaker=0
    while dif != 0 and breaker <= 40: #auto stop at 40 iterations
        r = 0
        for edge in G.edges:
            u,v,k = edge
            # The nodal pressure before correction
            bf = G.nodes[v].get('pressure',0)
            # Calculate a potential edge loss or use the given one
            edgeloss = calculate_edgeloss(G.edges[edge]['diameter'],G.edges[edge].get('c',G.edges[edge].get('roughness',c)),G.edges[edge]['length'],G.edges[edge].get("velocity",calculate_velocity(G.edges[edge]['diameter'],G.edges[edge].get('flowrate',G.edges[edge].get("flow",ValueError))))) # Calculate the edge loss
            # Calculate the new pressure
            G.nodes[v]['pressure'] = max(float(G.nodes[u]['pressure']) + float(G.nodes[u]['elevation']) - float(G.nodes[v]['elevation']) - float(G.edges[edge].get('edgeloss',edgeloss)) + G.edges[edge].get("pumphead",0), G.nodes[v]['pressure'])
            # THe corrected pressure
            af = G.nodes[v]['pressure']
            # Assign new head
            G.nodes[v]['head'] = af + G.nodes[v]['elevation']
            # The difference between actual and previous.
            r = max(r,af-bf)
        dif = r
        breaker+=1


def find_source_node(G):
    """get a list of all special nodes (Tanks, Reservoirs)

    Args:
        G (osmnx graph): The osmnx graph to search

    Returns:
        list: a list of nodes numbers with 'special' attributes
    """
    return list(nx.get_node_attributes(G,"special").keys())

def find_closest_source(G):
    """Find the closest source to every node

    Args:
        G (osmnx graph): The osmnx graph to check

    Returns:
        dict: A dictionary where keys are nodes and values are the distance to their closest source node.
    """
    sources = find_source_node(G)
    return nx.multi_source_dijkstra(G,sources)[0]
  

def make_boxplot(epa_G,char="pressure_dif"):
    """Make a boxplot showcasing a nodal char value based on distance from a source

    Args:
        epa_G (osmnx graph): The graph with nodal char
        char (str, optional): The nodal char to plot. Defaults to "pressure_dif".
    """
    # Find all the distances
    sources = find_source_node(epa_G)
    distances = nx.multi_source_dijkstra(epa_G,list(nx.get_node_attributes(epa_G,"special").keys()))[0]
    
    # Ignore source error
    for node in sources:
        del distances[node]
    errors = {}
    
    # Calculate the errors at all distances
    differences = nx.get_node_attributes(epa_G,char)
    for n,d in distances.items():
        if d not in errors:
            errors[d] = []    
        
        errors[d].append(differences[n])

    # Find the total error
    t_error = calculate_pressure_dif(epa_G,char)

    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(errors, orient='index').transpose()
    
    # Create boxplots
    fig, ax = plt.subplots(figsize=(30, 8))
    linewidth = dict(linewidth=2)
    df.boxplot(ax=ax,boxprops=linewidth,capprops=linewidth,whiskerprops=linewidth)
    ax.set_xlabel('Distance')
    ax.set_ylabel(char.capitalize() + ' Difference (m)')
    plt.title(t_error)
    plt.show()


def linear_interpolate(diameters, costs, diameter):
    """    Linearly interpolate to find the cost for a given diameter.
    Useful when comparing costs of network with non-standard pipe costs
    Or costs/diameters not in costs/diameters dictionary.


    Args:
        diameters (dict): A list of diameter sizes
        costs (dict): A list of costs of len(diameters)
        diameter (float): The diameter of interest

    Returns:
        float: A total cost
    """
    for i in range(len(diameters) - 1):
        if diameters[i] <= diameter <= diameters[i + 1]:
            # Linear interpolation formula
            d1, d2 = diameters[i], diameters[i + 1]
            c1, c2 = costs[i], costs[i + 1]
            interpolated_cost = c1 + (c2 - c1) * (diameter - d1) / (d2 - d1)
            return interpolated_cost
    # If diameter is outside the range, extrapolate
    if diameter < diameters[0]:
        # Extrapolate below the smallest known diameter
        d1, d2 = diameters[0], diameters[1]
        c1, c2 = costs[0], costs[1]
    else:
        # Extrapolate above the largest known diameter
        d1, d2 = diameters[-2], diameters[-1]
        c1, c2 = costs[-2], costs[-1]
    
    extrapolated_cost = c1 + (c2 - c1) * (diameter - d1) / (d2 - d1)
    return extrapolated_cost

ft_to_m = 0.3048

def calc_cost(G, costs,ft=False):
    """Calculate the total cost of a network based on pipes

    Args:
        G (osmnx network): The network to calculate the costs of
        costs (dict): a dictionary with costs and diameters {"costs":[cost1,cost2,costn],"diameters":[diam1,diam2,diamn]}
        ft (bool, optional): Whether the network is in ft. Defaults to False.

    Returns:
        float: A total network cost
    """
    diameters = costs["diameters"]
    cost_per_meter = costs["cost"]

    for edge in G.edges(data=True):
        u, v, data = edge
        length = data.get('length', 0)  # Assuming length is in meters
        diameter = data.get('diameter', 0)  # Assuming diameter is in meters
        if ft:
            length = length * ft_to_m
        # Find or interpolate the cost per meter for the given diameter
        if diameter in diameters:
            edge_cost_per_meter = cost_per_meter[diameters.index(diameter)]
        else:
            edge_cost_per_meter = linear_interpolate(diameters, cost_per_meter, diameter)

        # Calculate the total cost for this edge
        total_cost = length * edge_cost_per_meter
        data['cost'] = total_cost
    
    return int(sum(nx.get_edge_attributes(G,"cost").values()))

