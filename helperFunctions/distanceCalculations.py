#***************************************************************************************#
# Title:        distanceCalculations.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       Functions pertaining to the calculation of distances between
#       FLC and treatment values
#***************************************************************************************#
from scipy import stats
import numpy as np
from scipy.spatial import distance

#Will calculate the distance between each patient with a choice of different
#distance methods to apply. The methods currently available are:
#   P = Pearson's correlation
#   S = Spearman's correlation
#   K = Kendall Tau correlation
#   E = Euclidean Distance

#kendall tau does not return 0 for correlation between the same patient, but instead something
#very close to 0. Because of this, we have to set a minimum bound on the kendall tau correlation
#after which it will become 0.
kendallTauTolerance = .001 #how low we are going to let the kendal tau correlation get
def calculateFLCDistance(normalizedMatrix, whichDist):
    numRows = normalizedMatrix.shape[0]
    distanceMatrix = np.empty([numRows, numRows])
    for row in range(0, numRows): #goes through each patient
        for column in range(0, numRows): #compares to every other patient
            if(whichDist == "P"):
                correlation = 1 - stats.pearsonr(normalizedMatrix[row,1:].astype(float), normalizedMatrix[column,1:].astype(float))[0]
            if(whichDist == "S"):
                correlation = 1 - stats.spearmanr(normalizedMatrix[row,1:].astype(float), normalizedMatrix[column,1:].astype(float))[0]                
            if(whichDist == "K"):
                correlation = 1 - stats.kendalltau(normalizedMatrix[row,1:].astype(float), normalizedMatrix[column,1:].astype(float))[0]
                if(correlation < kendallTauTolerance):
                    correlation = 0                            
            if(correlation > 1 or np.isnan(correlation)):
                correlation = 1
            distanceMatrix[row,column] = correlation
            distanceMatrix[column, row] = correlation
    return distanceMatrix

            
        