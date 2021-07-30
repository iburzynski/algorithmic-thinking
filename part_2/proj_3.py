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
import alg_cluster


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
    For a list of clusters, checks the pairwise distance between all possible 
    cluster pairs. Returns a tuple with the minimum distance and the pair of 
    clusters in ascending order by index.
    Base case function for fast_closest_pair.
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
    Helper function for fast_closest pair. 
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
        (min_dist_r, *r_points) = fast_closest_pair(r_side)
        # add midpoint index to right points to restore their original indices
        for point in r_points:
            point += mid
        # set return values with the closest pair from both sides
        (min_dist, *points) = min(l_pair, (min_dist_r, *r_points))
        # calculate location of the center line
        strip = (clusters[mid-1].horiz_center() + clusters[mid].horiz_center())
        c_line = strip / 2
        # get closest pair within the center strip
        s_pair = closest_pair_strip(clusters, c_line, min_dist)
        # return the closest pair from all three subsets
        return min((min_dist, *points), s_pair)
 
    
######################################################################
# Code for hierarchical clustering


def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function may mutate cluster_list
    
    Input: List of clusters, integer number of clusters
    Output: List of clusters whose length is num_clusters
    """
    
    return []


######################################################################
# Code for k-means clustering

    
def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function may not mutate cluster_list
    
    Input: List of clusters, integers number of clusters and number of iterations
    Output: List of clusters whose length is num_clusters
    """

    # position initial clusters at the location of clusters with largest populations
            
    return []


### Testing ###
# Cluster(fips_codes, horiz_pos, vert_pos, population, risk)
test1 = closest_pair_strip([alg_cluster.Cluster(set([]), -4.0, 0.0, 1, 0), alg_cluster.Cluster(set([]), 0.0, -1.0, 1, 0), alg_cluster.Cluster(set([]), 0.0, 1.0, 1, 0), alg_cluster.Cluster(set([]), 4.0, 0.0, 1, 0)], 0.0, 4.123106)
test2 = closest_pair_strip([alg_cluster.Cluster(set([]), 0, 0, 1, 0), alg_cluster.Cluster(set([]), 1, 0, 1, 0), alg_cluster.Cluster(set([]), 2, 0, 1, 0), alg_cluster.Cluster(set([]), 3, 0, 1, 0)], 1.5, 1.0)
# expected: one of the tuples in set([(1.0, 1, 2)])
print(test1)
print(test2)