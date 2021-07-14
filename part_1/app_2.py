"""
Code for Application portion of Module 2
"""

# general imports
import requests
import random
import time
import math
import pandas as pd
import matplotlib.pyplot as plt

from proj_1 import make_complete_graph
from proj_2 import compute_resilience, compute_resilience_bfs
from app_1 import DataGraph, ERGraph, PAGraph, compute_neighbors, compute_probability
from upa_trial import UPATrial

############################################
# Provided code

def copy_graph(graph):
    """
    Make a copy of a graph
    """
    new_graph = {}
    for node in graph:
        new_graph[node] = set(graph[node])
    return new_graph

def delete_node(ugraph, node):
    """
    Delete a node from an undirected graph
    """
    neighbors = ugraph[node]
    ugraph.pop(node)
    for neighbor in neighbors:
        ugraph[neighbor].remove(node)
    
def targeted_order(ugraph):
    """
    Compute a targeted attack order consisting of nodes of maximal degree.
    Returns a list of nodes.
    """
    # copy the graph
    new_graph = copy_graph(ugraph)
    
    order = []    
    while len(new_graph) > 0:
        max_degree = -1
        for node in new_graph:
            if len(new_graph[node]) > max_degree:
                max_degree = len(new_graph[node])
                max_degree_node = node
        
        neighbors = new_graph[max_degree_node]
        new_graph.pop(max_degree_node)
        for neighbor in neighbors:
            new_graph[neighbor].remove(max_degree_node)

        order.append(max_degree_node)
    return order

##########################################################
# Code for loading computer network graph

NETWORK_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_rf7.txt"

##########################################################
# Question 1
##########################################################

def plot_resilience(graphs, targeted=False):
    # Simulate attacks and store resilience lists in dataframe
    data = pd.DataFrame()
    for graph in graphs:
        # print(graphs[graph].get_edges())
        data[graph] = graphs[graph].attack(targeted=targeted)
    attack_type = "Targeted" if targeted else "Random"
    title = f"Comparison of Graph Resilience for {attack_type} Order"
    upa_m = graphs['upa'].get_m()
    uer_p = graphs['uer'].get_prob()
    # Plot the resilience of each graph
    plt.figure()
    data.plot()
    plt.title(title)
    plt.xlabel("Number of Nodes Removed")
    plt.ylabel("Size of Largest Connected Component")
    plt.legend(["Computer Network", f"ER Graph (p={uer_p})", 
               f"UPA Graph (m={upa_m})"])
    plt.savefig(f"{attack_type.lower()}_resilience.png")

    return data

def is_resilient(data, attack_size=0.2, resilience_threshold=0.25):
    total_nodes = data.shape[0] - 1
    attack_threshold = math.ceil(total_nodes * attack_size)
    nodes_remaining = total_nodes - attack_threshold
    answers = {}
    for graph in data.columns:
        cc_size = data.at[attack_threshold, graph]
        answers[graph] = cc_size >= (1 - resilience_threshold) * nodes_remaining 

    return answers

# Create example graphs
graphs = {"network": DataGraph(NETWORK_URL, directed=False)}
num_nodes = graphs["network"].get_nodes()
prob = compute_probability(graphs["network"], rounding="ceil")
num_neighbors = compute_neighbors(graphs["network"], rounding="ceil")

graphs["uer"] = ERGraph(num_nodes, prob, directed=False)
graphs["upa"] = PAGraph(num_nodes, num_neighbors, directed=False)

random_df = plot_resilience(graphs, targeted=False)

print(is_resilient(random_df))

print(graphs["network"].fast_targeted_order())
# print(graphs["network"].get_degree(136))    
# print(graphs["network"].get_graph()[22])





