"""
Student template code for Project 3
Student will implement five functions:

slow_closest_pair(cluster_list)
fast_closest_pair(cluster_list)
closest_pair_strip(cluster_list, horiz_center, half_width)
hierarchical_clustering(cluster_list, num_clusters)
kmeans_clustering(cluster_list, num_clusters, num_iterations)

where cluster_list is a 2D list of clusters in the plane
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
        
    def merge_clusters(self, other_cluster):
        """
        Merges one cluster into another. The merge uses the relative 
        populations of each cluster in computing a new center and risk.
        
        Note: this method mutates self!
        """
        if len(other_cluster.fips_codes()) == 0:
        
            return self
        
        else:
            self._fips_codes.update(set(other_cluster.fips_codes()))
 
            # compute weights for averaging
            self_weight = float(self._total_population)                        
            other_weight = float(other_cluster.total_population())
            self._total_population = self._total_population + other_cluster.total_population()
            self_weight /= self._total_population
            other_weight /= self._total_population
                    
            # update center and risk using weights (refactor if possible)
            self._vert_center = self_weight * self._vert_center + other_weight * other_cluster.vert_center()
            self._horiz_center = self_weight * self._horiz_center + other_weight * other_cluster.horiz_center()
            self._averaged_risk = self_weight * self._averaged_risk + other_weight * other_cluster.averaged_risk()
        
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
        
        # compute error as weighted squared dist from counties to cluster center
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
    Helper function that computes Euclidean distance between two clusters in a list

    Input: cluster_list is list of clusters, idx1 and idx2 are integer indices for two clusters
    
    Output: tuple (dist, idx1, idx2) where dist is distance between
    cluster_list[idx1] and cluster_list[idx2]
    """
    return (cluster_list[idx1].distance(cluster_list[idx2]), min(idx1, idx2), max(idx1, idx2))


def slow_closest_pair(clusters):
    """
    Compute the distance between the closest pair of clusters in a list (slow)

    Input: cluster_list is the list of clusters
    
    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.       
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
    half_width is the half the width of the strip (i.e; the maximum horizontal distance
    that a cluster can lie from the center line)

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] lie in the strip and have minimum distance dist.       
    """
    # initialize default return values
    (min_dist, point1, point2) = (float('inf'), -1, -1)

    # list indices of clusters with x vals within half_w distance of center line  
    strip = [idx for idx, clstr in enumerate(clusters) if (abs(clstr.horiz_center() - h_center) < half_w)]
    # sort indices by vertical position (ascending)
    strip.sort(key = lambda cluster: clusters[cluster].vert_center())
    # get number of matching clusters to set iteration ranges
    num_clusters = len(strip)
    
    # for each cluster in the strip up to the second-to-last...
    for idx_u in range(0, num_clusters - 1):
        # check its pairwise distance with the next strip clusters (maximum 3)
        for idx_v in range(idx_u + 1, min(idx_u + 3, num_clusters - 1) + 1):
            dist, idx1, idx2 = pair_distance(clusters, strip[idx_u], strip[idx_v])
            # update return variables if distance is less than current minimum
            if dist < min_dist:
                (min_dist, point1, point2) = (dist, idx1, idx2)

    return min_dist, point1, point2

def fast_closest_pair(clusters):
    """
    Compute the distance between the closest pair of clusters in a list (fast)

    Input: cluster_list is list of clusters SORTED such that horizontal positions of their
    centers are in ascending order
    
    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.       
    """
    num_clusters = len(clusters)
    if num_clusters <= 3:
        return slow_closest_pair(clusters)
    
    else:
        # split cluster list in half
        mid = num_clusters // 2
        l_side = clusters[0:mid]
        r_side = clusters[mid:]
        # find closest pair on left side
        l_pair = fast_closest_pair(l_side)
        # find closest pair on right side
        #py3 (min_dist_r, *r_points) = fast_closest_pair(r_side)
        (min_dist_r, r_point1, r_point2) = fast_closest_pair(r_side) #py2
        # add midpoint index to right points to restore their original indices
        #py3 for point in r_points:
        #py3    point += mid
        r_point1 += mid #py2
        r_point2 += mid #py2
        # set return values with the closest pair from both sides
        #py3 (min_dist, *points) = min(l_pair, (min_dist_r, *r_points))
        (min_dist, point1, point2) = min(l_pair, (min_dist_r, r_point1, r_point2)) #py2
        # calculate location of the center line
        strip = (clusters[mid-1].horiz_center() + clusters[mid].horiz_center())
        c_line = strip / 2
        # get closest pair within the center strip
        s_pair = closest_pair_strip(clusters, c_line, min_dist)
        # return the closest pair from all three subsets
        #py3 return min((min_dist, *points), s_pair)
        return min((min_dist, point1, point2), s_pair) #py2

############################################################
# Code to create sequential clustering
# Create alphabetical clusters for county data

def sequential_clustering(singleton_list, num_clusters, **kwargs):
    """
    Take a data table and create a list of clusters
    by partitioning the table into clusters based on its ordering
    
    Note that method may return num_clusters or num_clusters + 1 final clusters
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


######################################################################
# Code for hierarchical clustering


def hierarchical_clustering(clusters, num_clusters, **kwargs):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function may mutate cluster_list
    
    Input: List of clusters, integer number of clusters
    Output: List of clusters whose length is num_clusters
    """
    while len(clusters) > num_clusters:
        clusters.sort(key = lambda cluster: cluster.horiz_center())
        (_dummy, pt1, pt2) = fast_closest_pair(clusters)
        clusters[pt1].merge_clusters(clusters[pt2])
        del clusters[pt2]

    return clusters


######################################################################
# Code for k-means clustering

    
def kmeans_clustering(clusters, num_clusters, num_iterations=5):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function may not mutate cluster_list
    
    Input: List of clusters, integers number of clusters and number of iterations
    Output: List of clusters whose length is num_clusters

    As you implement KMeansClustering, here are a several items to keep in mind. 
    In line 4, you should represent an empty cluster as a Cluster object whose 
    set of counties is empty and whose total population is zero. The cluster 
    centers μf, computed by lines 2 and 8-9, should stay fixed as lines 5-7 are 
    executed during one iteration of the outer loop. To avoid modifying these 
    values during execution of lines 5-7, you should consider storing these 
    cluster centers in a separate data structure. Line 7 should be implemented 
    using the merge_clusters method from the Cluster class. merge_clusters will 
    automatically update the cluster centers to their correct locations based on 
    the relative populations of the merged clusters.
    """

    # position initial centers at locations of clusters with largest populations
    sort_key = lambda cluster: cluster.total_population()
    old_clusters = sorted(clusters, key=sort_key, reverse=True)[:num_clusters]

    for _dummy_i in range(num_iterations):
        # initialize a set of k (num_clusters) empty clusters
        new_clusters = [Cluster(set(), 0, 0, 0, 0) for _cl in range(num_clusters)]
        for cls in clusters:
            # compute distances between the current cluster and each center 
            dists = [(cls.distance(cent), ix_f) for ix_f, cent in enumerate(old_clusters)]
            # identify the closest center and merge the current cluster with it
            _dist, closest = min(dists)
            new_clusters[closest].merge_clusters(cls)
        old_clusters = new_clusters

    return new_clusters