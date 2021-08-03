# import project algorithms and cluster classes
from proj_3 import fast_closest_pair, slow_closest_pair
from clusters import Cluster, Clustering 

# general imports
import random
import timeit
import copy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns


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
    # create algorithms dictionary to select function based on alg argument
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

def plot_distortion(datasets):
    """
    Takes a list of integer values (corresponding to the numbers of counties in 
    the unified cancer data spreadsheets) and produces hierarchical and 
    k-means clusterings with sizes ranging from 6 to 20 clusters. 
    Then calculates the distortion values for each cluster method at each 
    cluster size and creates a .png plot of the distortion for each dataset.
    """
    for dataset in datasets:
        cluster_sizes = list(range(6, 21, 1))
        # compute hierarchical clustering distortions
        h_clusts = []
        for idx, num_clusters in enumerate(reversed(cluster_sizes)):
            if len(h_clusts) == 0:
                # make initial clustering using the entire dataset
                h_clusts.append(Clustering(dataset, 'hierarchical', 
                                           num_clusters, mute=True))
            else:
                # then use prior clustering's list to prevent redundant computation
                cluster_list = copy.deepcopy(h_clusts[idx - 1].get_cluster_list())
                h_clusts.append(Clustering(dataset, 'hierarchical', 
                                           num_clusters, cluster_list, 
                                           mute=True))
        h_distort = [cluster.compute_distortion(mute=True) for cluster in reversed(h_clusts)]

        # compute k-means clustering distortions
        k_clusts = [Clustering(dataset, 'k-means', size, mute=True) for size in cluster_sizes]
        k_distort = [clust.compute_distortion(mute=True) for clust in k_clusts]
        zipped = list(zip(cluster_sizes, h_distort, k_distort))
        df = pd.DataFrame(zipped, columns=['num_clusters', 'hierarchical', 'k-means'])
        df.set_index('num_clusters', inplace=True)
        
        # plot the dataframe
        fig, ax = plt.subplots(1, 1)
        ax.plot(df.index, df['hierarchical'], label="Hierarchical Clustering")
        ax.plot(df.index, df['k-means'], label="K-means Clustering")
        # change offset formatter to scientific notation
        formatter = mticker.ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-3,2))
        ax.yaxis.set_major_formatter(formatter)
        # remove offset from y-axis
        ax.yaxis.offsetText.set_visible(False)
        # redraw the canvas
        fig.canvas.draw()
        # put offset text in y-axis label
        offset = ax.yaxis.get_major_formatter().get_offset()
        plt.ylabel(f"Distortion ({offset})")
        plt.xlabel("Number of Clusters")        
        plt.title(f"Distortion for Hierarchical and K-means Clustering ({dataset} Counties)")
        # add legend and reverse the legend order to match the plot positions
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], frameon=False)
        sns.despine()
        plt.savefig(f'plots/distortion_{dataset}.png')


##########################################################
# Question 1
##########################################################

# plot_tests(50)

##########################################################
# Question 2
##########################################################

# h_clust_3108 = Clustering(num_counties=3108, cluster_method='hierarchical', 
#                           num_clusters=15)
# h_clust_3108.run_viz()

##########################################################
# Question 3
##########################################################

# k_clust_3108 = Clustering(num_counties=3108, cluster_method='k-means', 
#                           num_clusters=15)
# k_clust_3108.run_viz()

##########################################################
# Question 5
##########################################################

# h_clust_111 = Clustering(num_counties=111, cluster_method='hierarchical', 
#                          num_clusters=9)
# h_clust_111.run_viz()

##########################################################
# Question 6
##########################################################

# k_clust_111 = Clustering(num_counties=111, cluster_method='k-means', 
#                          num_clusters=9)
# k_clust_111.run_viz()

##########################################################
# Question 7
##########################################################

# h_distortion = h_clust_111.compute_distortion()
# k_distortion = k_clust_111.compute_distortion()


##########################################################
# Question 10
##########################################################

plot_distortion([111, 290, 896])