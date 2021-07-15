"""
Code for Application portion of Module 1
"""
# import URL for generating citation graph
from constants_1 import CITATION_URL
# import graph classes
from graphs import CitationGraph, ERGraph, PAGraph
# import helper functions
from helpers import compute_neighbors, compute_probability

##########################################################
# Question 1
##########################################################

citation_graph = CitationGraph(CITATION_URL)
citation_graph.plot_dist()

##########################################################
# Question 2
##########################################################

# Create/plot a small sample graph with ER algorithm
er_graph_50 = ERGraph(50, .5)
er_graph_50.plot_dist()

# Create/plot a random graph with ER algorithm comparable to citation graph
num_nodes = citation_graph.get_nodes()
prob = compute_probability(citation_graph)
er_graph_full = ERGraph(num_nodes, prob, directed=True)
er_graph_full.plot_dist()

##########################################################
# Question 4
##########################################################

# Create/plot a DPA graph with similar n and m values to the citation graph
num_neighbors = compute_neighbors(citation_graph)
dpa_graph = PAGraph(num_nodes, num_neighbors, directed=True)
dpa_graph.plot_dist()