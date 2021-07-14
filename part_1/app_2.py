"""
Code for Application portion of Module 2
"""

# general imports
import requests
import random
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from proj_1 import make_complete_graph
from proj_2 import compute_resilience, compute_resilience_bfs
from app_1 import DataGraph, ERGraph, PAGraph, compute_neighbors, compute_probability
from upa_trial import UPATrial

############################################
# Provided code
    
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
    axes = plt.gca()
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
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

def time_test(num_neighbors=5, runs=10):
    def test_targeted(graph, fast=False):
        """
        Runs a single test on a graph using targeted_order. Returns time value.
        """
        to_start = time.time()
        if fast:
            graph.fast_targeted_order()
        else:
            graph.targeted_order()
        to_end = time.time()
        to_time = to_end - to_start 
        return to_time
    
    def test_series(graph, fast=False, runs=10):
        time_list = []
        for run in range(runs):
            time_list.append(test_targeted(graph, fast=fast))

        return np.mean(time_list)

    time_test = pd.DataFrame()
    num_nodes = []
    to_times = []
    fto_times = []

    for step in range(10, 1000, 10):
        num_nodes.append(step)
        graph = PAGraph(step, num_neighbors, directed=False)

        # Test targeted_order
        to_times.append(test_series(graph, fast=False, runs=runs))
        # Test fast_targeted_order
        fto_times.append(test_series(graph, fast=True, runs=runs))

    time_test["num_nodes"] = num_nodes
    time_test["targeted_order"] = to_times
    time_test["fast_targeted_order"] = fto_times
    time_test.set_index("num_nodes", inplace=True)

    plt.figure()
    time_test.plot()
    plt.title("Regular vs. Fast Computation of Targeted Order (Desktop Python)")
    axes = plt.gca()
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    plt.xlabel(f"Size (n) of UPA Graph (m = {num_neighbors})")
    plt.ylabel("Running Time (Seconds)")
    plt.savefig(f"targeted_comparison_{runs}.png")
    
    return time_test

##########################################################
# Question 1
##########################################################

# Create graphs
graphs = {"network": DataGraph(NETWORK_URL, directed=False)}
num_nodes = graphs["network"].get_nodes()
prob = compute_probability(graphs["network"], rounding="ceil")
num_neighbors = compute_neighbors(graphs["network"], rounding="ceil")

graphs["uer"] = ERGraph(num_nodes, prob, directed=False)
graphs["upa"] = PAGraph(num_nodes, num_neighbors, directed=False)

random_df = plot_resilience(graphs, targeted=False)

##########################################################
# Question 2
##########################################################

print(is_resilient(random_df))

##########################################################
# Question 3
##########################################################

# time_df = time_test(runs=150)

##########################################################
# Question 4
##########################################################

targeted_df = plot_resilience(graphs, targeted=True)

##########################################################
# Question 5
##########################################################

print(is_resilient(targeted_df))



