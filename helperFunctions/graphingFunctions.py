#***************************************************************************************#
# Title:        graphingFunctions.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       Functions for graphing clusters/subclusters.
#***************************************************************************************#
import csv
import matplotlib.pyplot as plt
import numpy as np

#Graphs the clusters in the folder ./clusters
def graphClusters(clusterMatrix, currentCluster, segLength):
    noProbabilities = np.delete(clusterMatrix, 2, 1) #gets rid of the probability column
    for row in range(0, noProbabilities.shape[0]):
        if(int(noProbabilities[row,1]) != -1):
            plt.plot(noProbabilities[row,2:].astype(float))
    plt.title("Cluster " + str(currentCluster))
    plt.ylabel("Normalized FLC value")
    plt.xlabel("28-day Segment")
    plt.savefig("./graphs/" + str(currentCluster) + ".png")
    plt.close()