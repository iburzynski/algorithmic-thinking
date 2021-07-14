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
# import helper classes for efficient versions of PA algorithms
from dpa_trial import DPATrial
from upa_trial import UPATrial
# import functions from Project 2
from proj_2 import compute_resilience, compute_resilience_bfs
from helpers import copy_graph, delete_node

CITATION_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"

###################################
# Code for Question 1

class Graph():
    def __init__(self, directed=True):
        """
        Makes an empty graph. Parent of CitationGraph, ERGraph and DPAGraph.
        """
        self._graph = dict()
        self._directed = directed
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
        Calculates the number of edges in the graph. Divides total by 2 if the
        graph is undirected.
        """
        graph_edges = {key: len(val) for key, val in self.get_graph().items()}
        total_edges = sum(graph_edges.values())
        if not self._directed:
            total_edges = total_edges // 2
        return total_edges

    def is_directed(self):
        return self._directed

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

    def random_order(self):
        """
        Returns a list of the graph's nodes in randomized order. Used for
        simulating attacks.
        """
        attack_list = list(self.get_graph().keys())
        random.shuffle(attack_list)
        return attack_list

    def targeted_order(self):
        """
        Compute a targeted attack order consisting of nodes of maximal degree.
        Returns a list of nodes.
        """
        # copy the graph
        new_graph = copy_graph(self.get_graph())
        
        order = []    
        while len(new_graph) > 0:
            max_degree = -1
            for node in new_graph:
                if len(new_graph[node]) > max_degree:
                    max_degree = len(new_graph[node])
                    max_degree_node = node
            
            # neighbors = new_graph[max_degree_node]
            # new_graph.pop(max_degree_node)
            # for neighbor in neighbors:
            #     new_graph[neighbor].remove(max_degree_node)

            delete_node(new_graph, max_degree_node)

            order.append(max_degree_node)
        return order

    def fast_targeted_order(self):
        """
        """
        num_nodes = self.get_nodes()
        degree_sets = {degree: set() for degree in range(num_nodes)}
        graph = copy_graph(self.get_graph())

        for node in graph:
            degree = len(graph[node])
            degree_sets[degree].add(node)
        attack_list = []

        for degree in range(num_nodes - 1, -1, -1):
            while len(degree_sets[degree]) > 0: 
                target = degree_sets[degree].pop()
                for neighbor in graph[target]:
                    nb_degree = len(graph[neighbor])
                    degree_sets[nb_degree].remove(neighbor)
                    degree_sets[nb_degree - 1].add(neighbor)
                attack_list.append(target)
                delete_node(graph, target)
        return attack_list

    def attack(self, targeted=False, algorithm="uf"):
        """
        Simulates an attack on the graph and returns a list of values 
        corresponding to the graph's resilience at each stage of the attack.
        """
        algs = {"uf": compute_resilience,
                "bfs": compute_resilience_bfs
               }
        graph = self.get_graph()
        if not targeted:
            attack_order = self.random_order()
        else:
            attack_order = self.fast_targeted_order()
        
        return algs[algorithm](graph, attack_order)


class DataGraph(Graph):
    def __init__(self, graph_url, directed=True):
        """
        Creates a graph given the URL for a text representation of it.
        """
        super().__init__(directed=directed)
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


class CitationGraph(DataGraph):
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
    def __init__(self, num_nodes, probability, directed=True):
        """
        Takes a number of nodes and probability and returns a random directed graph
        as a dictionary.
        """
        self._num_nodes = num_nodes
        self._probability = probability
        self._directed = directed
        nodes = list(range(0, num_nodes))
        graph = {node: set() for node in nodes}
        if directed:
            for node in nodes:
                for neighbor in nodes:
                    if neighbor != node and random.random() < probability:
                        graph[node].add(neighbor)
        else:
            for node in range(num_nodes - 1):
                for neighbor in range(node+1, num_nodes):
                    if random.random() < probability:
                        graph[node].add(neighbor)
                        graph[neighbor].add(node)

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

class PAGraph(Graph):
    def __init__(self, num_nodes, num_neighbors, directed=True):
        self._num_nodes = num_nodes
        self._num_neighbors = num_neighbors
        self._directed = directed
        if directed:
            self._graph_type = "dpa"
            self._trial_obj = DPATrial(num_neighbors) 
        else:
            self._graph_type = "upa"
            self._trial_obj = UPATrial(num_neighbors)
        # make a complete graph with m nodes
        graph = make_complete_graph(num_neighbors)

        for node in range(num_neighbors, num_nodes):
            possible_neighbors = self._trial_obj.run_trial(num_neighbors)
            graph[node] = possible_neighbors
            if not directed:
                # Add reciprocal edges if graph is undirected
                for neighbor in graph[node]:
                    graph[neighbor].add(node)
        self._graph = graph

    def get_m(self):
        """
        Helper function to retrieve the graph's m value.
        """
        return self._num_neighbors

    def get_type(self):
        """
        Helper function to retrieve the graph's type.
        """
        return self._graph_type

    def plot_dist(self):
        """
        Generates the PA distribution plot as a .png file.
        """
        file_name = f"{self.get_type()}_plot_{self.get_nodes()}"
        title = f"Log/Log Distribution of {self.get_type().upper()}"
        title += f" Graph (n={self.get_nodes()}, m={self.get_m()})"
        super().plot_dist(file_name, title)

        return None

###################################
# Code for Output

# Question 1
############
# citation_graph = CitationGraph(CITATION_URL)
# citation_graph.plot_dist()

# # Question 2
# ############
# # Create/plot a small sample graph with ER algorithm
# # er_graph_50 = ERGraph(50, .5)
# # er_graph_50.plot_dist()

# # Calculate appropriate values for n and probability based on the citation graph
def compute_neighbors(graph, rounding="ceil"):
    round_type = {"ceil": math.ceil,
                  "floor": math.floor
                 }
    num_neighbors = round_type[rounding](graph.get_edges() / graph.get_nodes())

    return num_neighbors

def compute_probability(graph, rounding="ceil", precision=4):
    num_nodes = graph.get_nodes()
    num_edges = graph.get_edges()
    multiplier = 1 if graph.is_directed() else 2
    prob = num_edges * multiplier / (num_nodes * (num_nodes - 1))
    return round(prob , precision)

# # Create/plot a random graph with ER algorithm comparable to citation graph
# num_nodes = citation_graph.get_nodes()
# prob = compute_probability(citation_graph)
# er_graph_full = ERGraph(num_nodes, prob, directed=True)
# er_graph_full.plot_dist()

# # Question 4
# ############

# # Create/plot a DPA graph with similar n and m values to the citation graph
# num_neighbors = compute_neighbors(citation_graph)
# dpa_graph = PAGraph(num_nodes, num_neighbors, directed=True)
# dpa_graph.plot_dist()