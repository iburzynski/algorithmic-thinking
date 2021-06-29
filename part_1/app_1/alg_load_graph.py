"""
Provided code for Application portion of Module 1

Imports physics citation graph 
"""

# general imports
import requests
from proj_1 import compute_in_degrees, in_degree_distribution
# Set timeout for CodeSkulptor if necessary
#import codeskulptor
#codeskulptor.set_timeout(20)


###################################
# Code for loading citation graph

CITATION_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"

def load_graph(graph_url):
    """
    Function that loads a graph given the URL
    for a text representation of the graph
    
    Returns a dictionary that models a graph
    """
    graph_file = requests.get(graph_url)
    graph_text = graph_file.text
    graph_lines = graph_text.split('\n')
    graph_lines = graph_lines[ : -1]
    
    print("Loaded graph with", len(graph_lines), "nodes")
    
    answer_graph = {}
    for line in graph_lines:
        neighbors = line.split(' ')
        node = int(neighbors[0])
        answer_graph[node] = set([])
        for neighbor in neighbors[1 : -1]:
            answer_graph[node].add(int(neighbor))

    return answer_graph


def normalize_dist(dist):
    """
    Takes an unnormalized distribution and returns a normalized distribution.
    """
    # get the number of nodes in the graph
    num_nodes = sum(dist.values())
    # create normalized distribution
    return {degree: count / num_nodes for degree, count in dist.items()}


# load the citation graph
citation_graph = load_graph(CITATION_URL)
# compute the in-degree distribution
unnorm_dist = in_degree_distribution(citation_graph)
# normalize the distribution
norm_dist = normalize_dist(unnorm_dist)

# print(norm_dist)


