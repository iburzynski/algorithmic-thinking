# Project/Application 1: Analysis of Citation Graphs
## Project Code (proj_1.py)
The Module 1 project file contains Python code that creates dictionaries corresponding to some simple examples of graphs, as well as two short functions that compute information about the distribution of the in-degrees for nodes in these graphs. These functions are utilized in the Application component of Module 1, where they are used to analyze the degree distribution of a citation graph for a collection of physics papers.
## Application Code (app_1.py)
The Module 1 application file combines mathematical analysis with the code from the project file to analyze the real-world problem of how scientific papers get cited. This application analyzes the structure of graphs generated by citation patterns from scientific papers to determine if cited papers are chosen randomly (from within the domain of the paper) or if there is some hidden pattern. It contrasts the in-degree distribution of a specific graph with those of graphs generated by two different random processes (ER and DPA algorithms).

### Notes on implementation:
* The solution uses an object-oriented approach, creating a generic graph class with child classes for generating and plotting the citation, ER, and DPA graphs. By storing the citation graph as an object, this implementation allows for information (i.e. nodes, edges) to be captured and subsequently used to generate comparable random graphs without hardcoding values for *n*, *p* or *m*. This allows for the application to be used with any similarly formatted citation dataset by simply modifying the URL constant. The graphs and plots produced by the ER and DPA algorithms will scale automatically based on the properties of the input data.
* The graph dictionaries are rendered as Pandas dataframes and plotted using Seaborn to produce more visually appealing plots.

### Question 1:
Compute the in-degree distribution for the citation graph, normalize it, and then create a log/log plot of the points in this normalized distribution.  

![Log/Log Citation Plot](./plots/citation_plot.png)

### Question 2:
Algorithm ER generates random undirected graphs. The following simple modification of the algorithm can be made to generate random directed graphs: for every ordered pair of distinct nodes (*i*,*j*), the modified algorithm adds the directed edge from *i* to *j* with probability *p*.

![ER Plot](./plots/er_plot_27770.png)

Consider the shape of the in-degree distribution for an ER graph and compare its shape to that of the physics citation graph:
* **Is the expected value of the in-degree the same for every node in an ER graph?** 
  Yes. The ER algorithm produces edges for each node with each other node using the same random process with a fixed probability of success p. This results in each node having the same expected value for its in-degree (EV = n*p).
* **What does the in-degree distribution for an ER graph look like?** 
  The in-degree distribution for ER graphs has a bump-shaped curve indicative of a binomial distribution. This curvature becomes more pronounced when the value for n is increased. Small values of n produce a more random distribution of points without a clear shape, while very large n values produce a very pronounced bump shape.
* **Does the shape of the in-degree distribution plot for ER look similar to the shape of the in-degree distribution for the citation graph?** 
  No. The in-degree distribution for graphs produced by the ER algorithm display a bump-shaped curve, which is characteristic of binomial distributions, while the in-degree distribution for the citation graph approximates a linear shape with a downward slope. This indicates that the process that produced the citation graph is not a purely random process like ER, since if it were the citation graph would have a binomial in-degree distribution similar in appearance to an ER graph.

### Question 3:
Next, consider a different process for generating synthetic directed graphs. In this process, referred to as *Algorithm DPA*, a random directed graph is generated iteratively, wherein each iteration a new node is created, added to the graph, and connected to a subset of the existing nodes. This subset is chosen based on the in-degrees of the existing nodes. To generate a random directed graph in this process, the user must specify two parameters: *n*, which is the final number of nodes, and *m* (where m???n), which is the number of existing nodes to which a new node is connected during each iteration. Notice that *m* is fixed throughout the procedure. 

The algorithm starts by creating a complete directed graph on *m* nodes, thengrows the graph by adding *n*-*m* nodes, where each new node is connected to *m* nodes randomly chosen from the set of existing nodes. As an existing node may be chosen more than once in an iteration, we eliminate duplicates (to avoid parallel edges); hence, the new node may be connected to fewer than *m* existing nodes upon its addition.

* **Choose values for *n* and *m* that yield a DPA graph whose number of nodes and edges is roughly the same to those of the citation graph.**  
n = 27770 (number of nodes in citation graph)  
m = 13 (computed by dividing number of edges in the citation graph by number of nodes and rounding up)

### Question 4:
Implement the DPA algorithm, compute a DPA graph using the values from Question 3, and then plot the in-degree distribution.
To improve the efficiency of the algorithm, the provided code contains a DPATrial class that maintains a list of node numbers and that contains multiple instances of the same node. If the number of instances of each node number is maintained in the same ratio as the desired probabilities, a call to random.choice() produces a random node number with the desired probability.

![DPA Plot](./plots/dpa_plot_27770.png)

### Question 5:
Compare the in-degree distribution for the citation graph to the in-degree distribution for the DPA graph as constructed in Question 4. In particular, consider whether the shape of the two distributions are similar and if so which of the following three phenomena might be the cause of the similarity:

1. The "six degrees of separation" phenomenon
2. The "rich gets richer" phenomenon
3. The "Hierarchical structure of networks" phenomenon

* **Is the plot of the in-degree distribution for the DPA graph similar to that of the citation graph?**  
Yes, the plot of the in-degree distribution for the DPA algorithm is highly similar to that of the citation graph. Both display a near-linear shape with a negative slope, as well as a fan-shaped scattering of data points in the lower right. This means in both graphs as the in-degree value increases (x-value), its probability (y-value) declines, and this linear relationship becomes less pronounced as values of in-degree increase.
* **Which one of the three social phenomena listed above mimics the behavior of the DPA process?** 
The "rich gets richer" social phenomenon mimics the behavior of the DPA process, and by extension the process behind the citation graph. These type of processes can be classified as "preferential attachment" processes, wherein the probability of the next outcome in a series receiving a particular value is proportional to the number of outcomes that already had that particular value. Preferential attachment processes are often discussed using terms like "cumulative advantage" and "popularity contests".\
In the DPA algorithm, connections are produced between new nodes and existing nodes based on the in-degrees of existing nodes. As opposed to the ER algorithm where each node has the same chance of receiving a connection from any other node, here existing nodes with higher in-degree values have a higher probability of success in receiving an edge from a new node (i.e. a citation from an existing paper). In the pseudocode for DPA this probability is expressed as (indeg(j) + 1)/(totindeg + |V|), where |V| is the current number of nodes. As the number of nodes increases iteratively, the overall probability of any existing node receiving an edge from a new node decreases (since both the total in-degree and number of nodes is increasing). However, existing nodes with higher in-degree values remain favored over peers with lower values. This exemplifies the "rich gets richer" phenomenon, as early nodes with larger amounts of incoming connections will be progressively favored and keep accumulating a relative advantage over later, less-connected nodes.
* **Could one of these phenomena explain the structure of the physics citation graph?** Yes, the "rich gets richer" phenomenon explains the structure of the citation graph as well. As a real-world example of the process behind the structure of the DPA graph, we can consider the way in which older papers which have been frequently cited would have higher visibility to writers of new papers than papers that haven't received as many citations, and thus be more likely to be cited again, leading to a cumulative advantage over time. 