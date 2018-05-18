#***************************************************************************************#
# Title:        labClustering.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Version History:
#       05/15/2018      First Release
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

#Will run a normalization on the data and then cluster. Users can choose
#which normalization to run with whichNorm. The current options for normalization
#are:
#   "M" = min-max normalization
#   "Z" = Z-score normalization
#   "U" = unit norm scaling
denominatorTolerance = 1e-6 # how low the max value can be
def getClustersOnTrend(toBeNormalized, whichDist, whichNorm):
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
    dc.calculateFLCDistance(toBeNormalized, whichDist)

#    deleterow = 0
#    distancenames = []
#    for i in range(0, toBeNormalized.shape[0]):
#        if (i < toBeNormalized.shape[0] - deleterow): 
#            distancenames.append(toBeNormalized[i][0])
#            if(toBeNormalized[i][0] == 'MM-326.1' or toBeNormalized[i][0] == 'MM-326.3' ):
#                   distancenames.remove(toBeNormalized[i][0])
#                   toBeNormalized = np.delete(toBeNormalized, i, 0)                   
#                   deleterow = deleterow + 1
#            
#            minFLC = min(toBeNormalized[i][1:].astype(float))
#            for j in range(1, len(toBeNormalized[i])):
#                toBeNormalized[i][j] = (toBeNormalized[i][j].astype(float) - minFLC.astype(float))
#            maxFLC = max(toBeNormalized[i][1:].astype(float))
#            
#            for j in range(1, len(toBeNormalized[i])):           
#                if(maxFLC.astype(float) != 0):
#                    toBeNormalized[i][j] = (toBeNormalized[i][j].astype(float))/maxFLC.astype(float)
#    justNorm = np.delete(toBeNormalized, 0, axis=1).astype(float)
#    pearCor = np.zeros((justNorm.shape[0], justNorm.shape[0]))
#    for i in range(0, justNorm.shape[0]):        
#        all_zeros = not np.any(justNorm[i,:])
#        if (all_zeros):
#            print(toBeNormalized[i,:])
#            print(justNorm[i,:])
#            print(i)
#            print(all_zeros)
#        for j in range(0, justNorm.shape[0]):
#            pearson = stats.pearsonr(justNorm[i,:], justNorm[j,:])
#            correlation = pearson[0]
#           pearCor[i, j] = correlation
#            pearCor[j, i] = correlation            
#   print(" is it finite??")
#    print(np.all(np.isfinite(pearCor)))
    #pearsonTrend = pearsonsCorr(justNorm)
#    distPearson = distanceMatrix(pearCor)
#    #clusters = initialHDBProcessing(toBeNormalized, "hdbscanBinnedData")
#    clusters = initialHDBProcessing(distPearson, "distanceMinMax", distancenames)
#    sortedClusters = clusters[clusters[:,1].astype(float).argsort()]
#    return sortedClusters