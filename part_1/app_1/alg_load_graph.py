"""
Code for Application portion of Module 1

Imports physics citation graph 
"""

# general imports
import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from proj_1 import compute_in_degrees, in_degree_distribution

###################################
# Code for Question 1

CITATION_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"

def load_graph(graph_url):
    """
    Function that loads a graph given the URL for a text representation of it.
    
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


def plot_dist(dist, file_name, plot_title, variables=["Degree", "Probability"]):
    """
    Takes a normalized distribution and plots it as a log-log plot.
    Generates a .png image.
    """
    sns.set()
    data = pd.DataFrame(dist.items(), 
                        columns=variables)
    sns_plot = sns.regplot(x=variables[0], y=variables[1], 
                           data=data,
                           scatter_kws={'alpha':0.33},
                           fit_reg=False)
    sns_plot.set(xscale="log", yscale="log", 
                 title=plot_title)
    plt.savefig(f"{file_name}.png")

    return None


def make_plot(graph, file_name, plot_title, variables):
    """
    """
    # compute the in-degree distribution
    unnorm_dist = in_degree_distribution(graph)
    # normalize the distribution
    norm_dist = normalize_dist(unnorm_dist)
    plot_dist(norm_dist, file_name, plot_title=plot_title, variables=variables)

    return None


def make_citation_plot():
    """
    Loads the citation data and generates the distribution chart as a .png file.
    """
    # load the citation graph
    citation_graph = load_graph(CITATION_URL)
    title = "Log/Log Distribution of High Energy Physics Paper Citations"
    vars = ['Number_of_Citations', 'Percentage_of_Papers']
    # make the plot
    make_plot(citation_graph, file_name="citation_plot", plot_title=title,
              variables=vars)

    return None

###################################
# Code for Question 2

def make_random_graph(num_nodes, probability):
    """
    Takes a number of nodes and probability and returns a random directed graph
    as a dictionary.
    """
    nodes = list(range(0, num_nodes))
    graph = {node: set() for node in nodes}
    for node in nodes:
        for neighbor in nodes:
            if neighbor != node and random.random() < probability:
                graph[node].add(neighbor)
    # print(graph)
    return graph


def make_random_plot(num_nodes, probability):
    """
    """
    # make a random directed graph
    random_graph = make_random_graph(num_nodes, probability)
    title = f"Log/Log Dist. of Random Directed Graph (n={num_nodes}, p={probability})"
    # plot the random graph
    make_plot(random_graph, f"random_plot_{num_nodes}", plot_title=title)

    return None

###################################
# Code for Output

# Question 1
############
make_citation_plot()

# Question 2
############
# make_random_plot(50, .5)
# make_random_plot(27000, .0005)