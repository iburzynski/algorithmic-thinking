# import project algorithms
from proj_3 import hierarchical_clustering, kmeans_clustering

# general imports
import math
import urllib
import matplotlib.pyplot as plt
from PIL import Image

###################################################
# Constants and Helper Functions

DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
COUNTIES = [111, 290, 896, 3108]
DATA_SOURCES = {num: f"{DIRECTORY}data_clustering/unifiedCancerData_{num}.csv" 
                        for num in COUNTIES}
MAP_URL = DIRECTORY + "data_clustering/USA_Counties.png"
# Define clustering methods
CLUSTERING_METHODS = {'hierarchical': hierarchical_clustering,
                      'k-means': kmeans_clustering}
# Define colors for clusters.  Display a max of 16 clusters.
COLORS = ['Aqua', 'Yellow', 'Blue', 'Fuchsia', 'Black', 'Green', 'Lime', 
          'Maroon', 'Navy', 'Olive', 'Orange', 'Purple', 'Red', 'Brown', 'Teal']

# Helper functions
def circle_area(pop):
    """
    Compute area of circle proportional to population
    """
    return math.pi * pop / (200.0 ** 2)


def plot_clusters(data_table, cluster_list, method, num_counties, 
                  draw_centers = False):
    """
    Creates a plot of clusters of counties
    """

    fips_to_line = {}
    for line_idx in range(len(data_table)):
        fips_to_line[data_table[line_idx][0]] = line_idx
     
    # load map image
    map_img = Image.open(urllib.request.urlopen(MAP_URL))

    # scale plot to get size similar to CodeSkulptor version
    xpixels, ypixels = map_img.size
    DPI = 60.0                  # adjust this constant to resize your plot
    xinch = xpixels / DPI
    yinch = ypixels / DPI
    plt.figure(figsize=(xinch,yinch))
    implot = plt.imshow(map_img)
   
    # draw the counties colored by cluster on the map
    if not draw_centers:
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.scatter(x = [line[1]], y = [line[2]], s =  circle_area(line[3]), lw = 1,
                            facecolors = cluster_color, edgecolors = cluster_color)

    # add cluster centers and lines from center to counties            
    else:
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.scatter(x = [line[1]], y = [line[2]], s =  circle_area(line[3]), lw = 1,
                            facecolors = cluster_color, edgecolors = cluster_color, zorder = 1)
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            cluster_center = (cluster.horiz_center(), cluster.vert_center())
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.plot( [cluster_center[0], line[1]],[cluster_center[1], line[2]], cluster_color, lw=1, zorder = 2)
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            cluster_center = (cluster.horiz_center(), cluster.vert_center())
            cluster_pop = cluster.total_population()
            plt.scatter(x = [cluster_center[0]], y = [cluster_center[1]], s =  circle_area(cluster_pop), lw = 2,
                        facecolors = "none", edgecolors = "black", zorder = 3)

    plt.title(f"{method.capitalize()} Clustering Map for {num_counties} Counties ({len(cluster_list)} clusters)", fontsize=24)
    # Remove axes
    plt.axis('off')

    plt.savefig(f'plots/{method}_clustering_map_{num_counties}.png')
