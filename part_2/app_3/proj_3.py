"""
Functions for Project component of Module 3
Author: Ian Burzynski
"""
import math

class Cluster:
    """
    Class for creating and merging clusters of counties.
    """
    def __init__(self, fips_codes, horiz_pos, vert_pos, population, risk):
        """
        Creates a cluster based on the model of a set of counties' data.
        """
        self._fips_codes = fips_codes
        self._horiz_center = horiz_pos
        self._vert_center = vert_pos
        self._total_population = population
        self._averaged_risk = risk
        
    def __repr__(self):
        """
        String representation, assuming the module is "clusters".
        """
        rep = "clusters.Cluster("
        rep += str(self._fips_codes) + ", "
        rep += str(self._horiz_center) + ", "
        rep += str(self._vert_center) + ", "
        rep += str(self._total_population) + ", "
        rep += str(self._averaged_risk) + ")"
        
        return rep

    def fips_codes(self):
        """
        Gets the cluster's set of FIPS codes.
        """
        return self._fips_codes
    
    def horiz_center(self):
        """
        Gets the averged horizontal center of cluster.
        """
        return self._horiz_center
    
    def vert_center(self):
        """
        Gets the averaged vertical center of the cluster.
        """
        return self._vert_center
    
    def total_population(self):
        """
        Gets the total population for the cluster.
        """
        return self._total_population
    
    def averaged_risk(self):
        """
        Gets the averaged risk for the cluster.
        """
        return self._averaged_risk
           
    def copy(self):
        """
        Returns a copy of a cluster.
        """
        copy_cluster = Cluster(set(self._fips_codes), self._horiz_center, 
                               self._vert_center, self._total_population, 
                               self._averaged_risk)
        
        return copy_cluster

    def distance(self, other_cluster):
        """
        Computes the Euclidean distance between two clusters.
        """
        vert_dist = self._vert_center - other_cluster.vert_center()
        horiz_dist = self._horiz_center - other_cluster.horiz_center()
        return math.sqrt(vert_dist ** 2 + horiz_dist ** 2)
        
    def merge_clusters(self, other):
        """
        Merges one cluster into another. The merge uses the relative 
        populations of each cluster in computing a new center and risk.
        
        Note: this method mutates self!
        """
        if len(other.fips_codes()) == 0:
            return self
        else:
            self._fips_codes.update(set(other.fips_codes()))
            # Compute weights for averaging
            s_weight = float(self._total_population)                        
            o_weight = float(other.total_population())
            self._total_population = int(s_weight + o_weight)
            s_weight /= self._total_population
            o_weight /= self._total_population     
            # Update center and risk using weights
            s_vc = self._vert_center
            s_hc = self._horiz_center
            s_ar = self._averaged_risk
            o_vc = other.vert_center()
            o_hc = other.horiz_center()
            o_ar = other.averaged_risk()

            self._vert_center = s_weight * s_vc + o_weight * o_vc
            self._horiz_center = s_weight * s_hc + o_weight * o_hc
            self._averaged_risk = s_weight * s_ar + o_weight * o_ar
        
            return self

    def cluster_error(self, data_table):
        """
        Input: data_table is the original table of cancer data used in creating 
        the cluster.
        
        Output: The error as the sum of the square of the distance from each 
        county in the cluster to the cluster center (weighted by its population)
        """
        # Build hash table to accelerate error computation
        fips_to_line = {}
        for line_idx in range(len(data_table)):
            line = data_table[line_idx]
            fips_to_line[line[0]] = line_idx
        # Compute error as weighted squared dist from counties to cluster center
        total_error = 0
        counties = self.fips_codes()
        for county in counties:
            line = data_table[fips_to_line[county]]
            singleton_cluster = Cluster(set([line[0]]), line[1], line[2], 
                                        line[3], line[4])
            singleton_dist = self.distance(singleton_cluster)
            error = (singleton_dist ** 2) * singleton_cluster.total_population()
            total_error += error
        
        return total_error

######################################################
# Code for closest pairs of clusters

def pair_distance(cluster_list, idx1, idx2):
    """
    Helper function that computes Euclidean distance between two clusters in a 
    list.
    Input: cluster_list is list of clusters, idx1 and idx2 are integer indices 
    for two clusters.
    Output: tuple (dist, idx1, idx2) where dist is distance between
    cluster_list[idx1] and cluster_list[idx2].
    """
    return (cluster_list[idx1].distance(cluster_list[idx2]), min(idx1, idx2), 
                                        max(idx1, idx2))


def slow_closest_pair(clusters):
    """
    Compute the distance between the closest pair of clusters in a list (slow).
    Input: cluster_list is the list of clusters.
    Output: tuple of the form (dist, idx1, idx2) where the centers of the 
    clusters cluster_list[idx1] and cluster_list[idx2] have minimum distance 
    dist.       
    """
    (min_dist, point1, point2) = (float('inf'), -1, -1)
    # check distances between all possible pairs to identify closest
    for idx_u, cluster_u in enumerate(clusters):
        for idx_v, cluster_v in enumerate(clusters):
            if idx_v != idx_u:
                 dist = cluster_u.distance(cluster_v)
                 if dist < min_dist:
                     min_dist, point1, point2 = (dist, idx_u, idx_v)

    return min_dist, point1, point2

def closest_pair_strip(clusters, h_center, half_w):
    """
    Helper function to compute the closest pair of clusters in a vertical strip
    Input: cluster_list is a list of clusters produced by fast_closest_pair
    horiz_center is the horizontal position of the strip's vertical center line
    half_width is the half the width of the strip (i.e; the maximum horizontal 
    distance that a cluster can lie from the center line)
    Output: tuple of the form (dist, idx1, idx2) where the centers of the 
    clusters cluster_list[idx1] and cluster_list[idx2] lie in the strip and have 
    minimum distance dist.       
    """
    def in_strip(cluster):
        """
        Helper function that takes a cluster as input and returns a boolean 
        value indicating whether its x val is within half_w distance of the 
        center line.
        """
        return abs(cluster.horiz_center() - h_center) < half_w

    # Initialize default return values
    (min_dist, point1, point2) = (float('inf'), -1, -1)
    # List indices of clusters with x vals within half_w distance of center line  
    strip = [idx for idx, cluster in enumerate(clusters) if in_strip(cluster)]
    # Sort indices by vertical position (ascending)
    strip.sort(key = lambda cluster: clusters[cluster].vert_center())
    # Get number of matching clusters to set iteration ranges
    num_clusters = len(strip)
    # For each cluster in the strip up to the second-to-last...
    for ix_u in range(0, num_clusters - 1):
        # Check its pairwise distance with the next strip clusters (maximum 3)
        for ix_v in range(ix_u + 1, min(ix_u + 3, num_clusters - 1) + 1):
            dst, idx1, idx2 = pair_distance(clusters, strip[ix_u], strip[ix_v])
            # Update return variables if distance is less than current minimum
            if dst < min_dist:
                (min_dist, point1, point2) = (dst, idx1, idx2)

    return min_dist, point1, point2

def fast_closest_pair(clusters):
    """
    Compute the distance between the closest pair of clusters in a list (fast)
    Input: cluster_list is list of clusters SORTED such that horizontal 
    positions of their centers are in ascending order.    
    Output: tuple of the form (dist, idx1, idx2) where the centers of the 
    clusters cluster_list[idx1] and cluster_list[idx2] have minimum distance 
    dist.       
    """
    num_clusters = len(clusters)

    if num_clusters <= 3:
        return slow_closest_pair(clusters)
    else:
        # Split the cluster list in half
        mid = num_clusters // 2
        l_side = clusters[0:mid]
        r_side = clusters[mid:]
        # Find the closest pair on the left side
        l_pair = fast_closest_pair(l_side)
        # Find the closest pair on the right side
        (min_dist_r, r_point1, r_point2) = fast_closest_pair(r_side) #py2
        # Add midpoint index to right points to restore their original indices
        r_point1 += mid
        r_point2 += mid
        # Set return values with the closest pair from both sides
        min_dist, point1, point2 = min(l_pair, (min_dist_r, r_point1, r_point2)) 
        # Calculate the location of the center line
        strip = (clusters[mid-1].horiz_center() + clusters[mid].horiz_center())
        c_line = strip / 2
        # Get closest pair within the center strip
        s_pair = closest_pair_strip(clusters, c_line, min_dist)
        # Return the closest pair from all three subsets
        return min((min_dist, point1, point2), s_pair)

def sequential_clustering(singleton_list, num_clusters, **kwargs):
    """
    Take a data table and create a list of clusters by partitioning the table 
    into clusters based on its ordering
    Note that method may return num_clusters or num_clusters + 1 final clusters.
    """
    cluster_list = []
    cluster_idx = 0
    total_clusters = len(singleton_list)
    cluster_size = float(total_clusters)  / num_clusters
    
    for cluster_idx in range(len(singleton_list)):
        new_cluster = singleton_list[cluster_idx]
        if math.floor(cluster_idx / cluster_size) != \
           math.floor((cluster_idx - 1) / cluster_size):
            cluster_list.append(new_cluster)
        else:
            cluster_list[-1] = cluster_list[-1].merge_clusters(new_cluster)
            
    return cluster_list

def hierarchical_clustering(clusters, num_clusters, **kwargs):
    """
    Compute a hierarchical clustering of a set of clusters.
    Note: the function may mutate cluster_list.
    Input: List of clusters, integer number of clusters.
    Output: List of clusters whose length is num_clusters.
    """
    while len(clusters) > num_clusters:
        clusters.sort(key = lambda cluster: cluster.horiz_center())
        (_dummy, pt1, pt2) = fast_closest_pair(clusters)
        clusters[pt1].merge_clusters(clusters[pt2])
        del clusters[pt2]

    return clusters

def kmeans_clustering(clusters, num_clusters, num_iterations=5):
    """
    Computes the k-means clustering of a set of clusters, returning a list of 
    clusters of length num_clusters.
    """
    # Position initial centers at locations of clusters with largest populations
    sort_key = lambda cluster: cluster.total_population()
    old = sorted(clusters, key=sort_key, reverse=True)[:num_clusters]

    for _idx in range(num_iterations):
        # Initialize a set of k (num_clusters) empty clusters
        new = [Cluster(set(), 0, 0, 0, 0) for _cluster in range(num_clusters)]
        for cls in clusters:
            # Compute distances between the current cluster and each center
            dists = [(cls.distance(cent), idx) for idx, cent in enumerate(old)]
            # Identify the closest center and merge the current cluster with it
            _dist, closest = min(dists)
            new[closest].merge_clusters(cls)
        old = new

    return new