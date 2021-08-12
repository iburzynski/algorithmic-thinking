"""
Code for Application component of Module 3
Author: Ian Burzynski
"""
# Import project algorithms and cluster classes
from proj_3 import fast_closest_pair, slow_closest_pair
from clusters import Cluster, Clustering 

# General imports
import copy
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import pandas as pd
import random
import seaborn as sns
import timeit

############################################
# Functions for answering Questions 1
############################################

def gen_random_clusters(num_clusters):
    """
    Creates a list of clusters where each cluster corresponds to one randomly 
    generated point in the square with corners (±1,±1).
    """
    cluster_list = []
    for cluster in range(num_clusters):
        x_pos = random.uniform(-1, 1)
        y_pos = random.uniform(-1, 1)
        cluster_list.append(Cluster(set(), x_pos, y_pos, 0, 0))

    return cluster_list

def closest_pair_time(size, iterations, alg='slow'):
    """
    Runs a specified number of timeit tests for a specified closest pair 
    algorithm and returns the average running time.
    """
    # Create algorithms dictionary to select function based on alg argument
    algs = {
        'slow': slow_closest_pair,
        'fast': fast_closest_pair
    }
    def closure():
        """
        Closure for running timeit test with provided arguments.
        """
        clusters = gen_random_clusters(size)
        return algs[alg](clusters)

    return timeit.timeit(closure, number = iterations) / iterations

def closest_pair_tests(iterations):
    """
    Runs timeit tests with slow and fast closest pair algorithms for cluster 
    lists of sizes 2 to 200 with the specified number of iterations.
    Returns the time results in a Pandas dataframe.
    """
    col = ['slow_closest_pair', 'fast_closest_pair']
    df = pd.DataFrame(columns=col)
    df.index.name = "Clusters"
    for size in range(2, 201):
        slow_time = closest_pair_time(size, iterations)
        fast_time = closest_pair_time(size, iterations, alg='fast')
        row = pd.Series([slow_time, fast_time], index=col, name=size)
        df = df.append(row)

    return df

def plot_tests(iterations):
    """
    Runs closest pair tests with the specified number of iterations, plots the 
    resulting dataframe and saves the plot as a .png file.
    """
    df = closest_pair_tests(iterations)
    plt.figure()
    df.plot()
    plt.title("Running Time of Slow vs. Fast Closest Pairs in Desktop Python")
    plt.xlabel("Number of Initial Clusters")
    plt.ylabel("Running Time in Seconds")
    plt.legend(frameon=False)
    sns.despine()
    plt.savefig('plots/closest_pair_tests.png')

    return None

############################################
# Function for answering Questions 10
############################################

def plot_distortion(datasets):
    """
    Takes a list of integer values (corresponding to the numbers of counties in 
    the unified cancer data spreadsheets) and produces hierarchical and 
    k-means clusterings with sizes ranging from 6 to 20 clusters. 
    Then calculates the distortion values for each cluster method at each 
    cluster size and creates a .png plot of the distortion for each dataset.
    """
    def make_k_clust(size):
        """
        Helper function to make a k-means clustering object of a given input 
        size. 
        """
        return Clustering(dataset, 'k-means', size, mute=True)

    for dataset in datasets:
        cluster_sizes = list(range(6, 21, 1))
        # Compute hierarchical clustering distortions
        h_clusts = []
        for idx, num_clusters in enumerate(reversed(cluster_sizes)):
            if len(h_clusts) == 0:
                # Make initial clustering using the entire dataset
                h_clusts.append(Clustering(dataset, 'hierarchical', 
                                           num_clusters, mute=True))
            else:
                # Use prior clustering's list to prevent redundant computation
                c_list = copy.deepcopy(h_clusts[idx - 1].get_cluster_list())
                h_clusts.append(Clustering(dataset, 'hierarchical', 
                                           num_clusters, c_list, mute=True))
        rev_hc = reversed(h_clusts)
        h_distort = [clust.compute_distortion(mute=True) for clust in rev_hc]
        # Compute k-means clustering distortions
        k_clusts = [make_k_clust(size) for size in cluster_sizes]
        k_distort = [clust.compute_distortion(mute=True) for clust in k_clusts]
        # Prepare the dataframe
        data = list(zip(cluster_sizes, h_distort, k_distort))
        cols = ['num_clusters', 'hierarchical', 'k-means']
        df = pd.DataFrame(data, columns=cols)
        df.set_index('num_clusters', inplace=True)
        # Plot the dataframe
        fig, ax = plt.subplots(1, 1)
        ax.plot(df.index, df['hierarchical'], label="Hierarchical Clustering")
        ax.plot(df.index, df['k-means'], label="K-means Clustering")
        # Change offset formatter to scientific notation
        formatter = mticker.ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-3,2))
        ax.yaxis.set_major_formatter(formatter)
        # Remove offset from y-axis
        ax.yaxis.offsetText.set_visible(False)
        # Redraw the canvas
        fig.canvas.draw()
        # Put offset text in y-axis label
        offset = ax.yaxis.get_major_formatter().get_offset()
        plt.ylabel(f"Distortion ({offset})")
        plt.xlabel("Number of Clusters")
        title = "Distortion for Hierarchical and K-means Clustering"
        title += f" ({dataset} Counties)"
        plt.title(title)
        # Add legend and reverse the legend order to match the plot positions
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], frameon=False)
        sns.despine()
        # Save the plot
        os.makedirs(f'./plots', exist_ok=True)
        plt.savefig(f'plots/distortion_{dataset}.png')

##########################################################
# Question 1
##########################################################

plot_tests(50)

##########################################################
# Question 2
##########################################################

h_clust_3108 = Clustering(num_counties=3108, cluster_method='hierarchical', 
                          num_clusters=15)
h_clust_3108.run_viz()

##########################################################
# Question 3
##########################################################

k_clust_3108 = Clustering(num_counties=3108, cluster_method='k-means', 
                          num_clusters=15)
k_clust_3108.run_viz()

##########################################################
# Question 5
##########################################################

h_clust_111 = Clustering(num_counties=111, cluster_method='hierarchical', 
                         num_clusters=9)
h_clust_111.run_viz()

##########################################################
# Question 6
##########################################################

k_clust_111 = Clustering(num_counties=111, cluster_method='k-means', 
                         num_clusters=9)
k_clust_111.run_viz()

##########################################################
# Question 7
##########################################################

h_distortion = h_clust_111.compute_distortion()
k_distortion = k_clust_111.compute_distortion()

##########################################################
# Question 10
##########################################################

plot_distortion([111, 290, 896])