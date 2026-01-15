import gurobipy as gp
from gurobipy import GRB
from math import pi
import math
import os

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

import osmnx as ox
import networkx as nx

class LinearProgram():
    """The Linear Program module.
    """
    def __init__(self, network, parameters, network_type):
        """This class can run different types of MILP dependign on need
        It can be run on different network types:
            - "full" : a full network with no cluster information
            - "cluster" : from the cluster algorithm with additional information regarding connections
            - "simple" : just like full but with no topology changes. All features set
        All of these can be run with set pump and water tower locations
        the head of these features can be optimized or set if already specified 

        To specify pump locations either:
            - give edges a set pumphead attribute
            - give edges a 'special': 'pump' in order to let program optimize the value
        
        Args:
            network (osmnx graph): the graph to optimize
            parameters (dict): the parameters needed to run the model, see examples
            network_type (str): What type of MILp to run, "full","cluster",or "simple"
        """
        
        # Load arguments
        self.G = network 
        self.parameters = parameters
        self.network_type = network_type

        # Check how much to print
        self.verbose = self.parameters.get('verbose')
        if self.verbose == None:
            self.verbose = 0
        
        # Total network demand must be 0 (sources - demand)
        net_dem = sum(list(nx.get_node_attributes(self.G, "demand").values()))
      
        if self.verbose:
            print("Demand Check")
            print("\tTotal demand:", round(net_dem,3))





        if round(net_dem,3) != 0:
            self.raise_error(ValueError,"total network demand must equal 0")
        
        # Prepare network
        # print("Preparing Network")
        self.prepare_network()
        
        # Showcase information
        if self.verbose:
            print("Network Statistics")
            print("\tNodes:",len(self.G)," ","Edges:",len(self.G.edges))

            print("Ready")

    def raise_error(self,error,desc):
        """an attempt at custom error messages

        Args:
            error (Error): the error type
            desc (str): the error description

        Raises:
            error: desc
        """
        raise error(desc)

    def prepare_network(self):
        """prepare the network (self.G) by ensuring all edges have equal and opposite edges, velocities, and that special nodes are designated as such.
        Also creates the gurobi model self.m
        """
        self.G = self.G.to_directed()
        for edge in list(self.G.edges): # Add opposite edges
            u,v,k = edge
            if "flowrate" in self.G.edges[edge]: #Q=VA, Q = V * D^2/4 * pi
                #import math
                self.G.edges[edge]['velocity'] = abs(self.G.edges[edge]["flowrate"] / 1000) / ((self.G.edges[edge]['diameter'] ** 2 / 4) * math.pi)

            if "velocity" not in self.G.edges[edge]:
                self.G.edges[edge]['velocity'] = self.parameters['velocity']
            if (v,u,k) not in self.G.edges:
                attrs = self.G.edges[edge].copy()
                if 'pumphead' in attrs:
                    attrs['pumphead'] *= -1
                self.G.add_edge(v,u,k,**attrs)

        for node in self.G.nodes: 
            if self.G.nodes[node]['demand'] < 0:
                spesh = self.G.nodes[node].get('special',0)
                if spesh in ["source",'Reservoir']:
                    self.G.nodes[node]['special'] = 'Reservoir'
                elif spesh in ['Tank','tower','tank']:
                    self.G.nodes[node]['special'] = "Tank"
        self.m = gp.Model()
        # print("\tPrepared")


    def load_vars(self):
        """Load the variables for the lp model
        """
        parameters = self.parameters


        # print("LOADING VARIABLES")
        ## Storage for the variables
        # For head loss
        self.edge_flow = {} #edge:edge_var, representing flow
        self.edge_loss = {} #edge:pressure_var, representing p loss over a pipe
        self.edge_loss_pos = {} # A positive component of edge_loss

        # For diameter sizes
        self.actual_diam = {} # edge: diam, The actual diameter of the edge
        self.edge_diam = {} # edge : {diamsize : binary}, Supports multiple sizes
        self.edge_diam_in = {} # node : Diam_size_var}, Will be constrained so var is >= all ins at that size (basically the size of the edge)

        # For equal and opposite edges
        self.edge_flow_up = {} # Posotive edge flow
        self.edge_flow_down = {} # negative edge flow
        self.up_down_flow = {} # binary edge_flow direction

        # For edge flows
        self.node_edges_var= {} #node: [edge_flow_vars], flows out

        # For node pressure
        self.node_pressure = {} #node:pressure_var, representing pressure at a node

        # For pump heads
        self.pump_head = {} # edge: pumphead, Pump head allowed, locations are not found, must be given in the Graph as pumphead

        ## Create the variables
        for node in self.G.nodes: # For every node
            n = self.m.addVar(name="NP"+str(node)) # Node Pressure

            # Node diameter 
            size = self.m.addVar(name="Dmaxin"+str(node),ub=max(parameters['diameters']))# Max size is the biggest diameter
            self.edge_diam_in[node] = size

            # All edges in and out
            edges = [] 
            

            ## Iterate throug edges attached to the node
            for edge in self.G.edges(node):
                u,v = edge

                # Add edge loss variables
                ef = self.m.addVar(name="EF"+str(edge[0])+"_"+str(edge[1]),lb=-float("inf"))
                el = self.m.addVar(name="EL"+str(edge[0])+"_"+str(edge[1]),lb=-float("inf"))
                elp = self.m.addVar(name="ELP"+str(edge[0])+"_"+str(edge[1]))

                # Add edge diameter options
                self.edge_diam[edge] = {} 
                for s in parameters['diameters']:
                    d = self.m.addVar(name="diam"+str(s)+"_"+str(edge[0])+"_"+str(edge[1]),vtype=GRB.BINARY) # make the variable for the size
                    self.edge_diam[edge][s] = d # Assign to the dict
                    
                ad = self.m.addVar(name="true_size"+str(edge[0])+"_"+str(edge[1]))
                self.actual_diam[edge] = ad
                
                # Incorporate pump head 
                special = self.G.edges[(edge[0],edge[1],0)].get('special')
                if special == "pump":
                    hd = self.G.edges[(u,v,0)]['pumphead']
                    if hd != None:
                        self.pump_head[edge] = float(hd)
                    else:
                        self.pump_head[edge] = self.m.addVar(name="pumphead"+str(edge[0])+"_"+str(edge[1]))
                else:
                    self.pump_head[edge] = 0

                # Edge flow
                efd = self.m.addVar(name="EFDown"+str(edge[0])+"_"+str(edge[1]))
                efu = self.m.addVar(name="EFUp"+str(edge[0])+"_"+str(edge[1]))
                edf = self.m.addVar(name="EDF"+str(edge[0])+"_"+str(edge[1]),vtype=GRB.BINARY)

                # Attach the variables to the storage for easy reference
                self.edge_flow[edge] = ef
                self.edge_loss[edge] = el
                self.edge_loss_pos[edge] = elp

                self.edge_flow_down[edge] = efd
                self.edge_flow_up[edge] = efu
                self.up_down_flow[edge] = edf

                edges.append(ef)

            self.node_edges_var[node] = edges
            self.node_pressure[node] = n
        self.m.update()

    def load_obj(self):
        """Create model objective. WIP for customization
        """
        parameters=self.parameters
        # print("LOADING OBJECTIVE")
        cost = gp.LinExpr()
       
        for edge in self.G.edges:
            u,v,k = edge
            for s in range(len(self.edge_diam[(u,v)])):
                cost += self.edge_diam[(u,v)][parameters['diameters'][s]]*parameters['cost'][s]*self.G.edges[edge]['length']
        # Other costs
            # cost += self.pump_head[(u,v)] * parameters["pump_cost_m"] 
        # for node in self.G.edges:
        #     cost+= self.tower_head[u] * parameters["tower_cost_m"]

        # Minimize the cost
        self.m.setObjective(
            cost,
            GRB.MINIMIZE
            )

        self.m.update()

    def hazen_williams(self,edge,length,velo,parameters,size):
        """The hazen-williams equation

        Args:
            length (float): length (m)
            velo (float): velocity (m/s)
            parameters (dict): the parameters used for roughness 'C'
            size (float): pipe diameter (m)

        Returns:
            float: the head loss of the edge
        """
        C = self.G.edges[edge].get('roughness',parameters.get('C',120))
        size_m = size / 1  # Convert diameter from mm to meters

        # head_loss = 10.67 * (length / (C**1.852)) * ((velo / size_m)**1.852)
        # return head_loss
        # parameters = self.parameters
        head_loss = 10.67 * (velo*(size_m/2)**2*pi)**1.852 * C**(-1.852) * (size_m)**-4.8704 * length
        return head_loss
        
    def load_edge_constr(self):
        """load all the base constraints for any LP
        """
        parameters = self.parameters

        # Edges we have already created constraints for
        seen = [] # For single edge constraints

        # Head Loss Edge Setup
        rijd = {} # The head loss
        ## Iterate through edges for flow
        for edge in self.G.edges:
            u,v,k = edge
            
            # Inverse Edge Flow
            if (u,v) not in seen or (v,u) not in seen: # Only need to do one or the other
                self.m.addConstr(self.edge_flow[(u,v)] == -1 * self.edge_flow[(v,u)],name="inverse"+str(u)+"_"+str(v))
                seen.append((u,v))
                seen.append((v,u))

            # Magnitude of Flow
            self.m.addConstr(self.edge_flow_up[(u,v)] - self.edge_flow_down[(u,v)] == self.edge_flow[(u,v)])

            # Available diameters 
            diam_both = gp.LinExpr() # Both up and down on an edge
            diam_up = gp.LinExpr() # up means u,v

            rijd[(u,v)] = {} # HEadloss

            # Calculate the true diameter size and the head loss
            diam_size = gp.LinExpr()
            tru_velo = gp.LinExpr() # Actual velo for the edge, must be less than a certain amt (min_velo)
            all_velo = []
            multiplier = gp.LinExpr()

            ## Iterate through all diameter options to find their head loss
            for size in parameters['diameters']:
                diam_both += self.edge_diam[(u,v)][size]
                diam_both += self.edge_diam[(v,u)][size]

                diam_up += self.edge_diam[(u,v)][size]

                # Set up headloss for each edge
                if "flowrate" in self.G.edges[edge]: #Q=VA, Q = V * D^2/4 * pi
                    velo = abs(self.G.edges[edge]["flowrate"] / 1000) / ((size**2 / 4) * pi)

                    # velo = abs(self.G.edges[edge]["flowrate"])/(self.G.edges[edge]['diameter']**2/4*pi)

                   
                else:
                    velo = float(self.G.edges[edge]['velocity'])

                all_velo.append(round(velo,5))
                tru_velo += round(velo,5)*self.edge_diam[(u,v)][size]
                

                # head_loss = 10.67 * (velo*(size/2)**2*pi)**1.852 * parameters['C']**(-1.852) * (size)**-4.8704 * self.G.edges[edge]['length']
                head_loss = self.hazen_williams(edge,self.G.edges[edge]['length'],velo,parameters,size)
                head_loss = round(head_loss,5)
                
                
                if edge[0] in nx.get_node_attributes(self.G,"special") or edge[1] in nx.get_node_attributes(self.G,'special'):
                    if self.pump_head[(u,v)] != 0: # If there is a pump head then the head loss is 0
                        head_loss = 0

                # all_h_loss[edge].append(head_loss)
                ####

                rijd[(u,v)][size] = gp.LinExpr(self.edge_diam[(u,v)][size]*head_loss) # Save the head loss
                
                diam_size += self.edge_diam[(u,v)][size]*size # Add to the diameter size
                multiplier += self.edge_diam[(u,v)][size] # for mutliplier, if all are 0 then diameter can be 0 (no pipe breaks the no size down constr)


            # print(all_velo,self.G.edges[edge]['flowrate'])
            # Ensure that the velocity is below the min_velo 
            if min(all_velo) >= parameters["min_velo"]: # Set velocity below a certain level
                self.m.addConstr(tru_velo<=min(all_velo)) # If the lowest velocity is still to high then just be lowest
            else: # Otherwise
                self.m.addConstr(tru_velo<=parameters['min_velo']) # If wiggle room then just be less than requirement

            self.m.addConstr(self.actual_diam[(u,v)] == diam_size)

            # If the entwork is the fire flow ensure we don't shrink pipes
            if self.network_type == "fire": # Fire flows
                if 'diameter' in self.G.edges[edge]: # Never shrink a pipe
                    self.m.addConstr(diam_size >= self.G.edges[edge]['diameter']*multiplier)
            
            self.m.addConstr(diam_both <= 1) # Only 1 type of pipe can exist on an edge

            # Max Flows
            self.m.addConstr(self.edge_flow_up[(u,v)]<=parameters["QMAX"],"Qmaxup"+str(u)+"_"+str(v))
            self.m.addConstr(self.edge_flow_down[(u,v)]<=diam_up*parameters["QMAX"],"Qmaxdown"+str(u)+"_"+str(v))

            # Min Flows
            self.m.addConstr(self.edge_flow_up[(u,v)]+self.edge_flow_down[(u,v)]>=diam_up*parameters["QMIN"],"QMin_up")

            # Constrain pos and neg flows so one is 0
            self.m.addConstr(self.edge_flow_down[(u,v)] <= (1-self.up_down_flow[(u,v)])*parameters["BigM_Q"]) #up down flow is an indicator of direction
            self.m.addConstr(self.edge_flow_up[(u,v)] <= self.up_down_flow[(u,v)]*parameters["BigM_Q"])
        ## iterate for head loss and pressures
        for edge in self.G.edges:
            u,v,k = edge
            # Pressure constraints
            self.m.addConstr(self.edge_loss_pos[(u,v)] == sum(rijd[(u,v)].values()))
           
            self.m.addConstr(self.edge_loss[(u,v)] == self.edge_loss_pos[(u,v)] - self.edge_loss_pos[(v,u)])

            # Bernoulis # Needed when pumps needed, otherwise fine
            self.m.addConstr(self.node_pressure[u] + self.G.nodes[u]['elevation'] + (1-(sum(self.edge_diam[(u,v)].values())))*parameters["BigM_P"] + self.pump_head[(u,v)] >=  self.node_pressure[v] + self.edge_loss[(u,v)] + self.G.nodes[v]['elevation'])
            # self.m.addConstr(self.node_pressure[u] + self.G.nodes[u]['elevation'] + self.tower_pressure[u] <=  self.node_pressure[v] + self.edge_loss[(u,v)] + self.G.nodes[v]['elevation'] + (1-(sum(self.edge_diam[(u,v)].values())))*parameters["BigM_P"] )
            
            # Max and min pressure constraints
            if self.G.nodes[u].get("special") not in ["Reservoir","Source","source"]: # If it isn't a source
                self.m.addConstr(self.node_pressure[u] >= parameters["MinP"])
                self.m.addConstr(self.node_pressure[u] <= parameters['MaxP'])
            else:
                self.m.addConstr(self.pump_head[(u,v)] + self.node_pressure[u] <= parameters["MaxP"]) # The source head and the pump head cannot be more than the max pressure
            
        
                
    def load_keep_edges_constr(self):
        """If all_edges in the parameters is set to True then all edges must have a diameter
        """
        for edge in self.G.edges:
            u,v,k = edge  
            if self.parameters.get("all_edges"): # If we are supposed to keep all edges
                self.m.addConstr(sum(self.edge_diam[(u,v)].values()) + sum(self.edge_diam[(v,u)].values()) == 1) # Force existing ONLY IF SPECIFIED


    def load_node_constr(self):
        """load all node based constraints, 
        """
        parameters = self.parameters
        # Nodal Constraints
        for node in self.G.nodes:
            #  Flow in - Flow out = Demand
            self.m.addConstr(float(self.G.nodes[node]['demand']) == sum(self.node_edges_var[node]),"flowrate"+str(node))

            
            if float(self.G.nodes[node]['demand']) >= 0: # If it has demand (excludes reservoirs and water towers)
                diam_in_each_size = [gp.LinExpr() for _ in range(len(parameters['diameters']))]
                for n in self.G.neighbors(node): #Everything in
                    for s in range(len(parameters['diameters'])):
                        diam_in_each_size[s] += self.edge_diam[(n,node)][parameters['diameters'][s]]
                
                for n in self.G.neighbors(node): 
                    for s in range(len(parameters['diameters'])):
                        self.m.addConstr(self.edge_diam[(node,n)][parameters['diameters'][s]] <= sum(diam_in_each_size[s:len(parameters['diameters'])]))
            else: # if the node is the source
                for edge in self.G.edges(node):
                    u,v = edge
                    self.m.addConstr(self.edge_diam[(u,v)][max(parameters['diameters'])] <= 1) 

            self.m.addConstr(self.node_pressure[node] + self.G.nodes[node]['elevation'] <= parameters['max_head'])
            

    # def load_set_max_head_constr(self):
    #     parameters = self.parameters
    #     # Nodal Constraints
        
        # i = 0
        # for node in nx.get_node_attributes(self.G,"special"):
        #     if i == 0:
        #         one = gp.LinExpr(self.node_pressure[node] + self.G.nodes[node]['elevation'])
        #     i = 1
        #     self.m.addConstr(one == self.node_pressure[node] + self.G.nodes[node]['elevation'])# == parameters['max_head'])

    def load_set_pressure_constr(self):
        """If a node has a set pressure that should not change use set_pressure in its attributes and call this
        """
        for node in self.G.nodes:
            if "set_pressure" in self.G.nodes[node]:
                self.m.addConstr(self.node_pressure[node] <= self.G.nodes[node]['set_pressure'])


    def load_cluster_constr(self):
        """When using the clustering algorithm there has to be extra constraints to ensure they are properly connected
        """
        #SPECIAL DIRECTIONALITY COSNTRAINTS
        
        for node in self.G.nodes:
            if "C" in str(node): # If the node is a cluster
                for edge in self.G.edges(node):
                    u,v = edge
                    if self.G.nodes[node]['demand'] <= 0: # it is a demand cluster
                        self.m.addConstr(sum(self.edge_diam[(v,u)].values()) == 0) # Force it to exist
                    else:
                        self.m.addConstr(sum(self.edge_diam[(u,v)].values()) == 0) # Force the opposite to exist
                        
        for edge in self.G.edges: 
            if 'exists' in self.G.edges[edge]:
                if "to" in self.G.edges[edge]['exists']:
                    direc = self.G.edges[edge]['exists'].split("to")
                    u = float(direc[0])
                    v = float(direc[1])
                    # self.m.addConstr(sum(self.edge_diam[(u,v)].values()) == 1)
                    self.m.addConstr(sum(self.edge_diam[(v,u)].values()) == 0)
                elif "both" == self.G.edges[edge]['exists']:
                    u,v,k = edge
                    self.m.addConstr(sum(self.edge_diam[(u,v)].values()) + sum(self.edge_diam[(v,u)].values()) == 1)
                    self.m.addConstr(self.edge_diam[(u,v)][max(self.parameters['diameters'])]+self.edge_diam[(v,u)][max(self.parameters['diameters'])]==1) 
        
        
        # Group edges to ensure clusters are connected and one is the right size
        groups = nx.get_edge_attributes(self.G,"group")
        g_sorted = {} 
        for edge,info in groups.items():
            u,v,k = edge
            if info not in g_sorted:
                g_sorted[info] = gp.LinExpr()
            g_sorted[info] += self.edge_diam[(u,v)][max(self.parameters['diameters'])] 
        if None in g_sorted:
            del g_sorted[None]
        for group in g_sorted.values():
            self.m.addConstr(group>=1) # One of these edges must be big


        # Force the edges to not be big (since predecessor edge was not)
        # print("ANTI_GROUP")
        anti_groups = nx.get_edge_attributes(self.G,"anti_group")
        g_sorted = {} 
        for edge,info in anti_groups.items():
            u,v,k = edge
            if info not in g_sorted:
                g_sorted[info] = gp.LinExpr()
            g_sorted[info] += self.edge_diam[(u,v)][max(self.parameters['diameters'])] 
        if None in g_sorted:
            del g_sorted[None]
        for group in g_sorted.values():
            self.m.addConstr(group==0) # NONE of these edges can be big

    def load_topo_constr(self):
        """Load the topography constraints from the parameters
        """
        # Topology Constraints
        # print("TOPOLOGY CONSTRAINTS") 
        e = gp.LinExpr()
        for edge in self.G.edges:
            u,v,k = edge
            e+=sum(self.edge_diam[(u,v)].values()) # ALl the edges

        n = len(self.G.nodes)

        # link density
        link_density = self.m.addVar(name="link_density")
        self.m.addConstr(link_density==2*e/(n*(n-1)))
        self.m.addConstr(self.parameters["link_density"]<=link_density)

        # node degree
        node_degree = self.m.addVar(name="node_degree")
        self.m.addConstr(node_degree==2*e/n)
        self.m.addConstr(self.parameters["node_degree"]<=node_degree)

        # meshedness
        meshedness = self.m.addVar(name="meshedness")
        self.m.addConstr(meshedness==(e-n+1)/(2*n-5))
        self.m.addConstr(self.parameters["meshedness"]<=meshedness)

    def load_source_out_constraints(self):
        """Sources can only push flow out.
        """

        for source in nx.get_node_attributes(self.G,'special'):
            for nbr in list(nx.neighbors(self.G,source)):
                self.m.addConstr(sum(self.edge_diam[(source,nbr)].values()) == 1)

    def load_constr(self):
        """load all the constraints based on network_type"""
        if self.network_type == "full":
            self.load_edge_constr()
            self.load_node_constr()
            self.load_topo_constr()
            self.load_set_pressure_constr()
           
        elif self.network_type == "cluster":
            self.load_edge_constr()
            self.load_node_constr()
            self.load_topo_constr()
            self.load_cluster_constr()

        elif self.network_type == "simple":
            self.load_edge_constr()
            self.load_node_constr()
            self.load_keep_edges_constr()
        
        elif self.network_type == "fire":
            self.load_edge_constr()
            self.load_node_constr()
            self.load_keep_edges_constr()

        else:
            self.raise_error(ValueError,"No Network Type Specified")


    def solve_model(self):
        """Solve the model
        """
        # print("SOVLING MODEL")
        
        self.m.update()
       
        # model settings
        self.m.setParam("LogToConsole", self.verbose) # Whether to show the dialogue
        self.m.setParam('MIPGap', self.parameters["Gap%"]) # Sets a limit on the difference between LP and MILP solutions
        self.m.setParam('TimeLimit', self.parameters['TimeLimit'])  # This sets a time limit
        # Set MIP focus to prioritize finding feasible integer solutions
        self.m.setParam('MIPFocus', self.parameters['MIPFocus'])

        # Solve
        self.m.optimize()

        # Save objective value
        self.obj_val = self.m.getObjective().getValue()

    def load_model(self):
        """Load the model constraints and variables"""
        self.load_vars()
        self.load_obj()
        self.load_constr()

    def run_model(self):
        """Run the model and save the output to the graph"""
        self.solve_model()
        self.assign_data()

    def get_objective_value(self):
        """Get the objective value"""
        return self.m.getObjective().getValue()

    def assign_data(self):
        """assign the model variable output to the edges and nodes
        Also remove edges with no diameter"""
        edges_to_remove = []
        
        for edge in self.G.edges:
            u,v,k = edge
            if round(sum(self.edge_diam[u,v].values()).getValue(),3) == 0:# round(self.actual_diam[(u,v)].x,6)
                edges_to_remove.append((u,v))
            else:
                self.G.edges[edge]['flowrate'] = round(self.edge_flow[(u,v)].x,6)
                self.G.edges[edge]['diameter'] = round(self.actual_diam[(u,v)].x,6)
                self.G.edges[edge]['edgeloss'] = round(self.edge_loss[(u,v)].x,6)
                self.G.nodes[u]['head'] = round(self.node_pressure[u].x+self.G.nodes[u]['elevation'],6)
                self.G.nodes[v]['head'] = round(self.node_pressure[v].x+self.G.nodes[v]['elevation'],6)
                self.G.edges[edge]['velocity'] = round(abs(self.G.edges[edge]["flowrate"]/1000)/(self.G.edges[edge]['diameter']**2/4*pi),6)

                self.G.nodes[u]['pressure'] = round(self.node_pressure[u].x,6)
                self.G.nodes[v]['pressure'] = round(self.node_pressure[v].x,6)
                self.G.nodes[u]['node_num'] = u
                self.G.nodes[v]['node_num'] = v
        
        self.G.remove_edges_from(edges_to_remove)

    def reload_G(self,G):
        """Update what the models graph should be

        Args:
            G (osmnx graph): the new graph
        """
        self.G = G

    def plot_network(self,G,edge_width=None,edge_color=None):
        """Plot the model network, 

        Args:
            G (osmnx graph): the graph to plot
            edge_width (str, optional): the parameter to graph on the edge. Defaults to None.
            edge_color (str, optional): the parameter to graph on the edge. Defaults to None.
        """
        nc = ["r" if G.nodes[node]['demand'] <= 0  else "k" for node in G.nodes] # Node colors are demand based
        ns = [abs(G.nodes[node]['demand']*400) if G.nodes[node]['demand'] >= 0 else 10 for node in G.nodes] 

        if edge_width:
            es = np.array(list(nx.get_edge_attributes(G,edge_width).values()))*10
        else:
            es = 1.5

        if edge_color:
            ecs = ox.plot.get_edge_colors_by_attr(G,edge_color,cmap='viridis_r')
        else:
            ecs = 'gray'

        ox.plot_graph(G,node_color=nc,node_size=ns,bgcolor=(0,0,0,0),edge_linewidth=es,edge_color=ecs,show=True)

    def super_plot(self,G,edge_label="velocity",node_label="name"):
        """Plot the graph but with labels

        Args:
            G (osmnx graph): the graph to plot
            edge_label (str, optional): the label for edges. Defaults to "velocity".
            node_label (str, optional): the node label. Defaults to "name".
        """
        fig,ax = plt.subplots(figsize=(25,20))

        x = nx.get_node_attributes(G,"x")
        y = nx.get_node_attributes(G,"y")
        pos = {}

        for n in x:
            pos[n] = [x[n],y[n]]
        node_labels = {node: int(G.nodes[node][node_label]) for node in G.nodes()} 
        nx.draw(
            G, pos, edge_color='black', width=1, linewidths=6,#[G2.edges[edge]['diam'] for edge in G2.edges],
            node_size=200, node_color='pink', alpha=0.9,
            labels=node_labels,arrowsize=20,arrows=True,font_size=20)

        elabels = {}
        for edge in G.edges:
            u,v,k = edge
            elabels[(u,v)] = round(G.edges[edge][edge_label],3)
          



        nx.draw_networkx_edge_labels(
            G, pos,
            edge_labels=elabels,
            font_color='red',
            font_size=20,
        )
        plt.show()
        
    