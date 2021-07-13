"""
Disjoint set class for find-union solution to compute_resilience in Project #2. 
"""
class DisjSet:
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