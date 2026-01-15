import osmnx as ox
import wntr
import shapely as sp
import networkx as nx
import pandas as pd
import numpy as np

from . import intersecting_method as im
from ..epanet_functions.iteration import * 



def wn_to_G(wn,crs="EPSG:4326"):
    """Convert a wntr network into an osmnx graph. 

    Args:
        wn (wntr network): a wntr network
        crs (str, optional): the projection for the model. Defaults to "EPSG:4326".

    Returns:
        osmnx graph: the preserved wntr network model
    """


   
    G = nx.MultiDiGraph()
    G.graph['crs'] = crs
    
    try:
        sim = wntr.sim.WNTRSimulator(wn)
        results = sim.run_sim()
        # print(results.node['demand'])
        
        fail = 0
    except:
        fail = 1

    if fail == 0:
        wn_gis = wntr.network.to_gis(wn,crs=crs)
        try:
            wn_gis.pipes["flowrate"] = abs(results.link['flowrate'].T)
        except:# Multiseries
            wn_gis.pipes['flowrate'] = abs(results.link['flowrate'].T[0])
        
        reservoir_inflow = results.node['demand'].loc[:, wn.reservoir_name_list] 
        tank_inflow = results.node['demand'].loc[:, wn.tank_name_list] 
    elif fail == 1:
        wn_gis = wntr.network.to_gis(wn,crs=crs)
        wn_gis.pipes["flowrate"] = 0
        reservoir_inflow = pd.DataFrame(zip(wn.reservoir_name_list,np.zeros(len(wn.reservoir_name_list)))).set_index(0).T
        tank_inflow = pd.DataFrame(zip(wn.tank_name_list,np.zeros(len(wn.tank_name_list)))).set_index(0).T
    # Add node
    
    for node_name, node in wn.nodes():
        attrs = node.to_dict()
        if isinstance(node, wntr.network.Junction):
            G.add_node(str(node_name), **attrs)
            dem = find_total_wntr_node_demand(wn.get_node(node))#attrs['demand_timeseries_list'][0]['base_val']
        elif isinstance(node, wntr.network.Reservoir):
            attrs['special'] = 'Reservoir'
            G.add_node(str(node_name), **attrs)
            dem = float(reservoir_inflow[str(node)].iloc[0])  #TODO: None is difficult for coded
        elif isinstance(node, wntr.network.Tank):
            attrs['special'] = 'Tank'
            G.add_node(str(node_name), **attrs)
            dem = float(tank_inflow[str(node)].iloc[0]) 
        G.nodes[str(node_name)]['x'] = attrs['coordinates'][0]
        G.nodes[str(node_name)]['y'] = attrs['coordinates'][1]
        G.nodes[str(node_name)]['demand'] = dem

    # Add edges
    for link_name, link in wn.links():
        attrs = link.to_dict()
        if isinstance(link, wntr.network.Pipe):
            G.add_edge(str(link.start_node.name), str(link.end_node.name), key=0,**attrs)
        elif isinstance(link, wntr.network.Pump):
            attrs['special'] = 'Pump'
            if 'pump_curve_name' in attrs and attrs['pump_curve_name'] in wn.curves:
                pump_curve = wn.get_curve(attrs['pump_curve_name'])
                # Set the pump head as the value at 0 flow (assuming a flat curve)
                if pump_curve.points:
                    attrs['pumphead'] = pump_curve.points[0][1]
            G.add_edge(str(link.start_node.name), str(link.end_node.name), key=0, **attrs)
        else:
            print("Edge Error")
            G.add_edge(str(link.start_node.name), str(link.end_node.name), key=0,**attrs)
        # elif isinstance(link, wntr.network.Valve): # TODO: remove?
        #     attrs['special'] = 'Valve'
        #     G.add_edge(int(link.start_node.name), int(link.end_node.name), key=0, **attrs)
        edge = (str(link.start_node.name), str(link.end_node.name), 0)
        s,e = tuple(im.calculate_geometry(G,edge[0],edge[1]).coords)
        
        geom = [s] + G.edges[edge]['vertices'] + [e]
        G.edges[edge]['geometry'] = sp.LineString(geom)
    nx.convert_node_labels_to_integers(G)
    for node in G.nodes:
        if 'demand_timeseries_list' in G.nodes[node]:
            del G.nodes[node]['demand_timeseries_list']
    return G
    
# def wntr_onto_osmnx(wn,raw_G):
#     """Given a wntr graph (wn) and networkx graph (G) that have the same topology place epanet simulation data onto graph G"""
#     #TODO: Include pumps and water towers
#     #TODO: Make it copy all features through a for loop
#     # TODO: Check what uses this and if I can stop
#     G = raw_G.copy()
#     sim = wntr.sim.WNTRSimulator(wn)
#     results = sim.run_sim()

#     wn_gis = wntr.network.to_gis(wn)
#     wn_gis.pipes["flowrate"] = abs(results.link['flowrate'].T)
#     wn_gis.junctions["pressure"] = abs(results.node['pressure'].T) 
#     wn_gis.junctions["demand"] = abs(results.node['demand'].T)

#     reservoir_inflow = results.node['demand'].loc[:, wn.reservoir_name_list] 
#     tank_inflow = results.node['demand'].loc[:, wn.tank_name_list] 

#     for node in G.nodes:
#         # if 'demand' not in G.nodes[node]:
#         #     G.nodes[node]['demand'] = 0.02
#         #     G.nodes[node]['pressure'] = 0
#         # else:
#         #     G.nodes[node]['demand'] = float(G.nodes[node]['demand'])
#         #     G.nodes[node]['pressure'] = float(G.nodes[node]['pressure'])

#         if str(node) in wn.node_name_list:
#             if wn.get_node(str(node)).node_type == "Junction":
#                 G.nodes[node]['demand'] = abs(float(wn_gis.junctions['demand'][str(node)]))
#                 G.nodes[node]['pressure'] = abs(float(wn_gis.junctions['pressure'][str(node)]))
#             if wn.get_node(str(node)).node_type == "Reservoir":
#                 dem = float(reservoir_inflow[str(node)].iloc[0]) 
#                 G.nodes[node]['demand'] = dem
#                 G.nodes[node]['special'] = wn.get_node(str(node)).node_type # TODO: Make nomenclature same throughout
#                 if dem >= 0:
#                     del G.nodes[node]['special']
#             elif wn.get_node(str(node)).node_type == "Tank":
#                 dem = float(tank_inflow[str(node)].iloc[0]) 
#                 G.nodes[node]['demand'] = dem
#                 G.nodes[node]['special'] = wn.get_node(str(node)).node_type # TODO: Make nomenclature same throughout
#                 if dem >= 0:
#                     del G.nodes[node]['special']

                
#     for edge in G.edges:
#         u,v,k= edge
        
#         k = str(u)+"_"+str(v)
#         k2 = str(v)+"_"+str(u)
#         if k in wn.link_name_list:
#             G.edges[edge]['flowrate'] = abs(float(wn_gis.pipes['flowrate'][k])) # Changed naming convention to flowrate
#         elif k2 in wn.link_name_list:
#             G.edges[edge]['flowrate'] = abs(float(wn_gis.pipes['flowrate'][k2]))
#         else:
#             G.edges[edge]['flowrate'] = 0
            
#     flowrates = results.link['flowrate']
    
#     for edge in G.edges:
#         u,v,k = edge
#         edge_types = ['pipe','pump']
#         for et in edge_types:
#             e1 = et+"_"+str(u) + "_" + str(v)+"_0" #TODO: Inconsistencies in how I name edges. pipe_u_v_k or just u_v
#             e2 = et+"_"+str(v) + "_" + str(u) +"_0"
            

#             if e1 in flowrates:
#                 G.edges[edge]['flowrate'] = float(flowrates[e1].iloc[0])
#             elif e2 in flowrates:
#                 G.edges[edge]['flowrate'] = float(flowrates[e2].iloc[0])

#     return G

    

# def osmnx_onto_wntr(G,wn):
#     """Transfer data from a G network to a wn network. Primarilty used for diameter transfer.


#     Args:
#         G (osmnx graph): taking data from this
#         wn (wntr network): and transfer to this

#     Returns:
#         wntr network: a modified wntr network
#     """
#     for edge in G.edges:
#         u,v,k = edge
#         edge_types = ['pipe','pump']
#         for et in edge_types:
#             e1 = et+"_"+str(u) + "_" + str(v)+"_0" #TODO: Inconsistencies in how I name edges. pipe_u_v_k or just u_v
#             e2 = et+"_"+str(v) + "_" + str(u) +"_0"
            
#             if e1 in wn.link_name_list:
#                 eoi = wn.get_link(e1)
#                 eoi.diameter = G.edges[edge]['diameter']
#             elif e2 in wn.link_name_list:
#                 eoi = wn.get_link(e2)
#                 eoi.diameter = G.edges[edge]['diameter']
            
#     return wn

def wn_to_inp(wn,filename):
    """Convert a wn network to a .inp file

    Args:
        wn (wntr network): network to save
        filename (str): location to save
    """
    wntr.network.io.write_inpfile(wn, filename)

def G_to_inp(G,filename,geometry=True):
    """Convert a G network to a .inp file 

    Args:
        G (osmnx graph): the osmnx graph to save
        filename (str): where to save network
        geometry (bool, optional): whether to keep geometry. Defaults to True.
    """
    wn_to_inp(G_to_wn(G,geometry),filename) # TODO: Need to work on how to define things as pumps

def wn_to_osm(wn,filename):
    """Convert a wn network to a .osm file

    Args:
        wn (wntr network): the wntr network to save
        filename (str): where to save the network
    """
    G_to_osm(wn_to_G(wn),filename)

def G_to_osm(G,filename):
    """Convert a G network to a .osm file

    Args:
        G (osmnx graph): the graph to save
        filename (str): where to save the graph
    """
    ox.save_graphml(G,filename)

def inp_to_wn(filename):
    """Convert a .inp file to a wn network

    Args:
        filename (str): the file location to load

    Returns:
        wntr network: the wntr network model stored in the inp file
    """
    wn = wntr.network.WaterNetworkModel(filename)    
    return wn

def inp_to_G(filename):
    """Convert a .inp file to a G network

    Args:
        filename (str): the file location to load

    Returns:
        osmnx graph: the osmnx graph converted from the wntr network inp file
    """
    wn = inp_to_wn(filename)
    G = wn_to_G(wn)
    return G

def osm_to_wn(filename):
    """ Convert a .osm file to a wn network

    Args:
        filename (filename): the file location to load

    Returns:
        wntr network: the wntr network converted from the osmnx graph osm file
    """
    return G_to_wn(osm_to_G(filename))

def osm_to_G(filename,edge_dtypes={"velocity":float,"flowrate":float,"headloss":float},node_dtypes={"demand":float,"head":float,"pressure":float,"street_count":float}):
    """Convert a .osm file to a G network. Reskin of  ox.load_graphml(). Note that when errors occur with loading try no default types and png.correct_attr_types()

    Args:
        filename (str): the file location of the graph
        edge_dtypes (dict, optional): what default edge types to use. Defaults to {"velocity":float,"flowrate":float,"headloss":float}.
        node_dtypes (dict, optional): what default node types to use. Defaults to {"demand":float,"head":float,"pressure":float}.

    Returns:
        _type_: _description_
    """
    G = ox.load_graphml(filename,edge_dtypes=edge_dtypes,node_dtypes=node_dtypes)
    correct_attr_type(G)
    return G


def load_real():

    # Load real Network
    real = inp_to_G("Input/Lakewood_no_valves.inp")
    wn_real = inp_to_wn("Input/Lakewood_no_valves.inp")
    # wn = png.G_to_wn(G2)
    real = transfer_flowrates(wn_real,real)





    def give_pumps_diams(G):
        for edge in G.edges:
            if not G.edges[edge].get('diameter'):
                G.edges[edge]['diameter'] = 0.6
        return G


    real = give_pumps_diams(real)
    # Make sure I convert to LPS and m
    for edge in real.edges: # Correct velocity  #TODO: Do this somewhere else
        real.edges[edge]['velocity'] = calculate_velocity(real.edges[edge]['diameter']*1000,real.edges[edge]['flowrate']/1000)

    for node in real.nodes:
        if 'elevation' in real.nodes[node]:
            real.nodes[node]['head'] = real.nodes[node].get('elevation') + real.nodes[node]['pressure']
    real.graph['crs'] = 'EPSG:2229'
    G_projected = ox.projection.project_graph(real, to_crs='EPSG:4326')  # NAD83 / California zone 6
    for node in G_projected.nodes:
        G_projected.nodes[node]['coordinates'] = (G_projected.nodes[node]['x'], G_projected.nodes[node]['y'])
    return G_projected
