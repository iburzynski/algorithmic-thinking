"""
Classes for creating and visualizing clusters of county-based cancer risk data 
and associated helper code.
Author: Ian Burzynski
"""
# Import helper code and cluster class
from helpers import CLUSTERING_METHODS, DATA_SOURCES, plot_clusters
from proj_3 import Cluster

# General imports
import requests

class Clustering:
    """
    Class for creating, plotting and comparing lists of clusters.
    """
    def __init__(self, num_counties, cluster_method, num_clusters, 
                 cluster_list=None, mute=False):
        """
        Initializes a clustering object, loads the data table corresponding 
        to the specified number of counties, and creates a cluster list using 
        the specified clustering method. For hierarchical clustering, an 
        existing cluster can be provided to create the cluster list with reduced 
        computation. The mute argument turns print statements on or off.
        """
        self._num_counties = num_counties
        self._cluster_method = cluster_method
        self._cluster_func = CLUSTERING_METHODS[cluster_method]
        self._num_clusters = num_clusters
        self._data_table = self.load_data_table(mute=mute)
        self._cluster_list = self.make_cluster_list(cluster_list)

    def load_data_table(self, mute=False):
        """
        Imports a table of county-based cancer risk datafrom a csv format file.
        """
        data_file = requests.get(DATA_SOURCES[self._num_counties])
        data = data_file.text
        data_lines = data.split('\n')

        if not mute:
            print(f"Loaded {len(data_lines)} data points")

        data_tokens = [line.split(',') for line in data_lines]
        return [[tokens[0], float(tokens[1]), float(tokens[2]), int(tokens[3]), 
                 float(tokens[4])] for tokens in data_tokens]

    def make_cluster_list(self, cluster_list):
        """
        Creates a list of clusters using the data table and clustering method 
        specified at initialization. For hierarchical clustering, if an existing 
        cluster list was provided at initialization, the existing list is used 
        to create the cluster list instead of beginning from the number of 
        counties in the data table. This feature reduces computation load when 
        creating multiple clusterings for comparison, such as in 
        plot_distortion().
        """
        if self._cluster_method == 'hierarchical' and cluster_list:
            return self._cluster_func(cluster_list, self._num_clusters)

        singleton_list = []
        for line in self._data_table:
            singleton_list.append(Cluster(set([line[0]]), line[1], 
                                  line[2], line[3], line[4]))
        
        return self._cluster_func(singleton_list, self._num_clusters)	

    def get_cluster_list(self):
        """
        Getter function for _cluster_list property.
        """
        return self._cluster_list

    def run_viz(self, draw_centers=True):
        """
        Plots the clustering on the US map graphic and saves as a .png file.
        """
        print(f"Plotting {self._num_clusters} {self._cluster_method} clusters")
        plot_clusters(self._data_table, self._cluster_list, 
                      self._cluster_method, self._num_counties, draw_centers)

        return None

    def compute_distortion(self, mute=False):
        """
        Computes the clustering's amount of distortion relative to its county 
        data table. 
        Mute argument turns print statement on or off. 
        """
        distortion = sum([clus.cluster_error(self._data_table) 
                          for clus in self._cluster_list])
        if not mute:
            method_str = self._cluster_method.capitalize()
            print(f"{method_str} Clustering Distortion: {distortion:.3E}")

        return distortion





  
        






        




