"""
Code for Application portion of Module 1
"""

# general imports
import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import math
# import functions from Project 1
from proj_1 import compute_in_degrees, in_degree_distribution, make_complete_graph
# import helper class for efficient version of DPA algorithm (Question 4)
from alg_dpa_trial import DPATrial

CITATION_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"

###################################
# Code for Question 1

class Graph():
    def __init__(self):
        """
        Makes an empty graph. Parent of CitationGraph, ERGraph and DPAGraph.
        """
        self._graph = dict()
        self._num_nodes = 0

        return None
    
    def get_graph(self):
        """
        Helper method for retrieving citation graph dictionary.
        """
        return self._graph

    def get_nodes(self):
        return self._num_nodes

    def get_edges(self):
        """
        Calculates the number of edges in the graph.
        """
        graph_edges = {key: len(val) for key, val in self.get_graph().items()}
        return sum(graph_edges.values())

    def get_raw_dist(self):
        """
        Calculates graph's raw in-degree distribution.
        """
        return in_degree_distribution(self.get_graph())

    def normalize_dist(self):
        """
        Creates normalized distribution from the graph's raw distribution.
        """
        return {deg: count / self._num_nodes for deg, count in self.get_raw_dist().items()}

    def plot_dist(self, file_name, plot_title, variables=["Degree", "Probability"]):
        """
        Creates a normalized distribution and plots it as a log-log plot.
        Generates a .png image.
        """
        sns.set()
        data = pd.DataFrame(self.normalize_dist().items(), 
                            columns=variables)
        sns_plot = sns.regplot(x=variables[0], y=variables[1], 
                            data=data,
                            scatter_kws={'alpha':0.33},
                            fit_reg=False)
        sns_plot.set(xscale="log", yscale="log", 
                    title=plot_title)
        plt.savefig(f"{file_name}.png")
        # Clear figure.
        plt.clf()
        return None


class CitationGraph(Graph):
    def __init__(self, graph_url):
        """
        Creates a citation graph given the URL for a text representation of it.
        """
        super().__init__()
        graph_file = requests.get(graph_url)
        graph_text = graph_file.text
        graph_lines = graph_text.split('\n')
        graph_lines = graph_lines[ : -1]
        
        # determine the number of nodes
        self._num_nodes = len(graph_lines)

        for line in graph_lines:
            neighbors = line.split(' ')
            node = int(neighbors[0])
            self._graph[node] = set([])
            for neighbor in neighbors[1 : -1]:
                self._graph[node].add(int(neighbor))
        
        self._num_edges = self.get_edges()
        print(f"Loaded graph with {self._num_nodes} nodes and {self._num_edges} edges")

        # create normalized distribution
        self._norm_dist = super().normalize_dist()

    def plot_dist(self):
        """
        Loads the citation data and generates the distribution chart as a .png file.
        """
        file_name = "citation_plot"
        title = "Log/Log Distribution of High Energy Physics Paper Citations"
        vars = ['Number_of_Citations', 'Percentage_of_Papers']
        super().plot_dist(file_name, title, vars)

        return None

###################################
# Code for Question 2

class ERGraph(Graph):
    def __init__(self, num_nodes, probability):
        """
        Takes a number of nodes and probability and returns a random directed graph
        as a dictionary.
        """
        self._num_nodes = num_nodes
        self._probability = probability
        nodes = list(range(0, num_nodes))
        graph = {node: set() for node in nodes}
        for node in nodes:
            for neighbor in nodes:
                if neighbor != node and random.random() < probability:
                    graph[node].add(neighbor)
        self._graph = graph

    def get_prob(self):
        """
        Helper method to get probability value used for edge creation.
        """
        return self._probability

    def plot_dist(self):
        """
        Generates the ER distribution plot as a .png file.
        """
        file_name = f"er_plot_{self.get_nodes()}"
        title = f"Log/Log Distribution of ER Graph (n={self.get_nodes()}, p={self.get_prob()})"
        super().plot_dist(file_name, title)

        return None

###################################
# Code for Question 4

class DPAGraph(Graph):
    def __init__(self, num_nodes, num_neighbors):
        self._num_nodes = num_nodes
        self._num_neighbors = num_neighbors
        # make a complete graph with m nodes
        complete_graph = make_complete_graph(num_neighbors)
        trial_obj = DPATrial(num_neighbors)

        for node in range(num_neighbors, num_nodes):
            complete_graph[node] = trial_obj.run_trial(num_neighbors)
        
        self._graph = complete_graph

    def get_m(self):
        """
        Helper function to retrieve the graph's m value.
        """
        return self._num_neighbors

    def plot_dist(self):
        """
        Generates the DPA distribution plot as a .png file.
        """
        file_name = f"dpa_plot_{self.get_nodes()}"
        title = f"Log/Log Distribution of DPA Graph (n={self.get_nodes()}, m={self.get_m()})"
        super().plot_dist(file_name, title)

        return None

###################################
# Code for Output

# Question 1
############
citation_graph = CitationGraph(CITATION_URL)
citation_graph.plot_dist()

# Question 2
############
# Create/plot a small sample graph with ER algorithm
# er_graph_50 = ERGraph(50, .5)
# er_graph_50.plot_dist()

# Calculate appropriate values for n and probability based on the citation graph
num_nodes = citation_graph.get_nodes()
num_edges = citation_graph.get_edges()
num_neighbors = math.ceil(num_edges / num_nodes)
prob = round(num_neighbors / num_nodes, 4)
# Create/plot a random graph with ER algorithm comparable to citation graph
er_graph_full = ERGraph(num_nodes, prob)
er_graph_full.plot_dist()

# Question 4
############

# Create/plot a DPA graph with similar n and m values to the citation graph
# dpa_graph = DPAGraph(num_nodes, num_neighbors)
# dpa_graph.plot_dist()