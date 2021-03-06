"""
Project functions for Application #1
"""
## import sample graphs for testing
# import constants_1

def make_complete_graph(num_nodes):
    """
    Takes a number of nodes and returns a dictionary corresponding to a 
    complete directed graph with the specified number of nodes.
    """
    graph = {}
    for node in range(num_nodes):
        graph[node] = set()
        for edge in range(num_nodes):
            if edge != node:
                graph[node].add(edge)
    return graph

def compute_in_degrees(digraph):
    """
    Takes a directed graph (represented as a dictionary) and computes the
    in-degrees for the nodes in the graph. Returns a dictionary with the
    same set of keys as digraph whose values are the number of edges whose
    head matches a particular node.
    """
    in_degrees = {}
    for node in digraph:
        if node not in in_degrees:
            in_degrees[node] = 0
        for edge in digraph[node]:
            if edge in in_degrees:
                in_degrees[edge] += 1
            else:
                in_degrees[edge] = 1

    return in_degrees

def in_degree_distribution(digraph):
    """
    Takes a directed graph (represented as a dictionary) and computes the
    unnormalized distribution of the in-degrees of the graph. Returns a 
    dictionary whose keys correspond to the in-degrees of nodes in the graph.
    """
    in_deg_dist = {}
    in_degrees = compute_in_degrees(digraph)

    for node in in_degrees:
        deg = in_degrees[node]
        if deg in in_deg_dist:
            in_deg_dist[deg] += 1
        else:
            in_deg_dist[deg] = 1

    return in_deg_dist