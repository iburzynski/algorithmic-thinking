"""
Project functions for Application #2
"""
from collections import deque

## import sample graphs for testing
# import constants_2 

# Code for challenge question
class DisjSet:
    """
    Disjoint set class for find-union solution to compute_resilience in 
    Project #2. 
    """
    def __init__(self, ugraph):
        """
        Initializes a disjointed set with no edges with nodes from an 
        undirected graph adjacency list.
        """
        self.ranks = {node: 0 for node in ugraph.keys()} # depth of the tree
        self.parents = {node: node for node in ugraph.keys()}
        # Create dict to record number of children (used to find largest cc)
        self.children = {node: 0 for node in ugraph.keys()}
  
    def find(self, x):
        """
        Takes a node and returns the representative of the set it is an element 
        of.
        """
        if (self.parents[x] != x):
            
            # if x is not the parent of itself
            # Then x is not the representative of
            # its set,
            self.parents[x] = self.find(self.parents[x])
            
            # so we recursively call Find on its parent
            # and move i's node directly under the
            # representative of this set  
        return self.parents[x] 
  
    def union(self, node_x, node_y):
        """
        Creates a union of two sets represented by node_x and node_y. 
        """     
        # Find current sets of x and y
        xset = self.find(node_x)
        yset = self.find(node_y)
  
        # Exit method if the nodes are already in the same set.
        if xset == yset:
            return
  
        # If ranks are different, put the lower ranked item under the higher.
        if self.ranks[xset] < self.ranks[yset]:
            self.parents[xset] = yset
            self.children[yset] += (1 + self.children[xset])
        else:
            self.parents[yset] = xset
            self.children[xset] += (1 + self.children[yset])
            # If ranks are same, move y under x and increment rank of x's tree.
            if self.ranks[xset] == self.ranks[yset]:
                self.ranks[xset] = self.ranks[xset] + 1

# Breadth-first search
def bfs_visited(ugraph, start_node):
    """
    Takes the undirected graph ugraph and the node start_node and returns the 
    set consisting of all nodes that are visited by a breadth-first search that 
    starts at start_node.
    """
    visited = set()
    # Initialize an empty double-ended queue object
    queue = deque()
    visited.add(start_node)
    queue.append(start_node)

    while len(queue) > 0:
        current_node = queue.pop()
        for neighbor in ugraph[current_node]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    
    return visited

# Connected components
def cc_visited(ugraph):
    """
    Takes the undirected graph ugraph and returns a list of sets, where each set
    consists of all the nodes (and nothing else) in a connected component, and 
    there is exactly one set in the list for each connected component in ugraph 
    and nothing else.
    """
    cc_list = []
    # Initialize a list of remaining nodes containing all nodes in the graph
    remaining = list(ugraph.keys())

    while len(remaining) > 0:
        # Select an arbitrary start node from the list of remaining nodes
        start_node = remaining[0]
        # Get all visited nodes from BFS beginning with that node
        visited = bfs_visited(ugraph, start_node)
        # Add the set of visited nodes to the CC list
        cc_list.append(visited)
        # Remove the visited nodes from the list of remaining nodes
        remaining = [item for item in remaining if item not in visited]

    return cc_list

def largest_cc_size(ugraph):
    """
    Takes the undirected graph ugraph and returns the size (an integer) of the 
    largest connected component in ugraph.
    """
    largest = 0
    cc_list = cc_visited(ugraph)

    for component in cc_list:
        size = len(component)
        if size > largest:
            largest = size

    return largest

# Graph resilience
def compute_resilience_bfs(ugraph, attack_order):
    """
    Takes an undirected graph and an ordered list of nodes to attack. For each 
    node in the list, the function removes the node and its edges from the graph 
    and then computes the size of the largest connected component for the 
    resulting graph. 
    Returns a list whose k + 1th entry is the size of the largest
    connected component in the graph after of the removal of the first k nodes 
    in attack_order. The first entry (indexed by zero) is the size of the 
    largest connected component in the original graph. 
    """
    resilience = [largest_cc_size(ugraph)]
    for attacked_node in attack_order:
        if attacked_node in ugraph:
            del ugraph[attacked_node]
            for edges in ugraph.values():
                edges.discard(attacked_node)
        resilience.append(largest_cc_size(ugraph))

    return resilience

# Challenge version using disjoint sets
def compute_resilience(ugraph, attack_order):
    """
    Efficient alternative to compute_resilience_bfs. Takes an undirected graph 
    ugraph, a list of nodes attack_order and iterates through the nodes in 
    attack order in reverse. For each node, the function reconstructs the ugraph
    as a disjoint graph based on its neighbors in the adjacency list.

    """
    def get_largest(disj):
        """
        Takes a disjoint graph and returns the size of the largest connected
        component (highest number of children + 1).
        """
        return max(disj.children.values()) + 1
    # Initialize the connected component size list as a deque to append leftwise
    cc_sizes = deque()
    # Create a disjoint set with nodes of ugraph
    disj = DisjSet(ugraph)

    # Build the post-attack graph with any nodes that weren't attacked
    safe_nodes = [node for node in ugraph if node not in attack_order]
    built_nodes = []
    for node in safe_nodes:
        for neighbor in ugraph[node]:
            if neighbor in built_nodes:
                disj.union(node, neighbor)
        built_nodes.append(node)
    if len(safe_nodes) == 0:
        # If every node is attacked, add 0 to end of sizes list
        cc_sizes.append(0)
    else:
        # If safe nodes exist, add size of largest connected component to list
        cc_sizes.append(get_largest(disj))
    
    # Connect the attacked nodes to their neighbors in reverse order
    for node in reversed(attack_order):
        for neighbor in ugraph[node]:
            if neighbor in built_nodes:
                disj.union(node, neighbor)
        # Get size of largest connected component by finding most common parent
        cc_sizes.appendleft(get_largest(disj))
        built_nodes.append(node)
    return list(cc_sizes)