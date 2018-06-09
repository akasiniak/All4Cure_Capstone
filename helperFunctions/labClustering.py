#***************************************************************************************#
# Title:        labClustering.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       All functions necessary for clustering of lab values. Clusters all patient
#       bins first by trend between consecutive lab values, and then further sub
#       clusters by scale. If the graph option is selected, the clusters/subclusters
#       will be graphed.
#***************************************************************************************#
from helperFunctions import distanceCalculations as dc
from scipy import stats
from sklearn import preprocessing
import numpy as np
import time
from hdbscan import HDBSCAN
import os
import glob
from helperFunctions import graphingFunctions as gr
#Will run a normalization on the data and then cluster. Users can choose
#which normalization to run with whichNorm. The current options for normalization
#are:
#   "M" = min-max normalization
#   "Z" = Z-score normalization
#   "U" = unit norm scaling
denominatorTolerance = 1e-6 # how low the max value can be
def getClustersOnTrend(toBeNormalized, whichDist, whichNorm, segLength, graph):
    if(whichNorm == "M"):
        for row in range(0, toBeNormalized.shape[0]): #go through every row in the matrix
            minFLC = min(toBeNormalized[row][1:].astype(float))
            for column in range(1, len(toBeNormalized[row])): #subtract the min from each index in the row
                toBeNormalized[row][column] = toBeNormalized[row][column].astype(float) - minFLC
            maxFLC = max(toBeNormalized[row][1:].astype(float))
            for column in range(1, len(toBeNormalized[row])): #divide the max from each index
                if(maxFLC.astype(float) > denominatorTolerance):
                    toBeNormalized[row][column] = (toBeNormalized[row][column].astype(float))/maxFLC
                else:
                    toBeNormalized[row][column] = denominatorTolerance
    if(whichNorm == "Z"):
        zscoreMatrix = stats.zscore(toBeNormalized[:,1:].astype(float), axis = 1)
        for row in range(0, toBeNormalized.shape[0]):
            for column in range(1, len(toBeNormalized[row])): #replace each index in toBeNormalized with Z-score
                if(np.isnan(zscoreMatrix[row][column - 1])):
                    toBeNormalized[row][column] = denominatorTolerance
                else:
                    toBeNormalized[row][column] = zscoreMatrix[row][column - 1]
    if(whichNorm == "U"):
        unitNorm = preprocessing.normalize(toBeNormalized[:,1:], norm = 'l2')
        for row in range(0, toBeNormalized.shape[0]):
            for column in range(1, len(toBeNormalized[row])): #replace each index in toBeNormalized with unit norm value
                toBeNormalized[row][column] = unitNorm[row][column - 1]
    distanceMatrix = dc.calculateFLCDistance(toBeNormalized, whichDist)
    clusters = clusterLabValues(distanceMatrix, toBeNormalized[:,0])

    #This is for creating the clustering files
    path = glob.glob('./clusters/*') #remove the cluster files currently in the folder
    for file in path:
        os.remove(file)
    sortedClusters = clusters[clusters[:,1].astype(float).argsort()]
    currentCluster = int(sortedClusters[0 , 1]) #gets the first available cluster
    clusterMatrix = np.empty(segLength + 3) #will contain the cluster data for each cluster
    for row in range(0, sortedClusters.shape[0]):
        if(currentCluster != int(sortedClusters[row, 1]) or row == sortedClusters.shape[0] - 1): #the current data point is in another cluster
            if(row == sortedClusters.shape[0] - 1): #for the last patient segment
                buildRow = np.empty(segLength + 3).astype(str)
                buildRow[0] = sortedClusters[row, 0]
                buildRow[1] = sortedClusters[row, 1]
                buildRow[2] = sortedClusters[row, 2]
                for row in range(0, toBeNormalized.shape[0]):
                    if(toBeNormalized[row, 0] == buildRow[0]):
                        for dataPoint in range(0, segLength):
                            buildRow[dataPoint + 3] = toBeNormalized[row,dataPoint + 1]
                clusterMatrix = np.vstack((clusterMatrix, buildRow))
            clusterMatrix = np.delete(clusterMatrix, 0, 0) #numpy.empty makes the first row random data so we get rid of it
            if(graph == 1):
                gr.graphClusters(clusterMatrix, currentCluster, segLength)
            np.savetxt("./clusters/cluster" + str(currentCluster) + ".csv", clusterMatrix.astype(str), delimiter = ',', fmt='%s')
            clusterMatrix = np.empty(segLength + 3)
            currentCluster = currentCluster + 1
        buildRow = np.empty(segLength + 3).astype(str)
        buildRow[0] = sortedClusters[row, 0]
        buildRow[1] = sortedClusters[row, 1]
        buildRow[2] = sortedClusters[row, 2]
        for row in range(0, toBeNormalized.shape[0]):
            if(toBeNormalized[row, 0] == buildRow[0]):
                for dataPoint in range(0, segLength):
                    buildRow[dataPoint + 3] = toBeNormalized[row,dataPoint + 1]
        clusterMatrix = np.vstack((clusterMatrix, buildRow))

#Clusters values based on FLC values. Takes in a distance matrix and a vecotr with all the corresponding
#patient numbers and outputs a cluster matrix
def clusterLabValues(workingMatrix, distanceNames):
    hdb_t1 = time.time()
    hdb = HDBSCAN(min_cluster_size=3, metric='precomputed').fit(workingMatrix)
    hdb_labels = hdb.labels_
    hdb_prob = hdb.probabilities_
    hdb_elapsed_time = time.time() - hdb_t1
    n_clusters_hdb_ = len(set(hdb_labels)) - (1 if -1 in hdb_labels else 0)
    merged = np.array(list(zip(distanceNames, hdb_labels, hdb_prob)))
    return(merged)