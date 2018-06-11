#***************************************************************************************#
# Title:        treatmentClustering.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       All functions necessary for clustering of treatment values.
#***************************************************************************************#
import numpy as np
import globalVariables as gv
from scipy import stats
from scipy.spatial import distance
from hdbscan import HDBSCAN
from sklearn.cluster import DBSCAN
from sklearn import metrics

# Accepts different option selections as parameters, and uses previously made med sequence matrix
# to create and return treatement distance matrix
def getDistancesFromMeds(noBits, signBitOption, aggregateTreatmentsOption, binaryDistanceMethod, useManualDistanceMethod, segmentLength):
    # Pull in previously made med sequence matrix and delete the patient identifiers
    treatmentData = np.delete(gv.medSequenceMatrix, 0, 1)
    # Create an array of the binary strings which are encoded with each segment's treatments
    monthlyBinaryTreatmentVectors = []
    for i in range(0, treatmentData.shape[0]):
        temp = []
        for j in range(0, segmentLength):
            binaryTreatment = convTreatmentToBinaryArray(treatmentData[i, j], signBitOption, noBits)
            temp.append(binaryTreatment)
        monthlyBinaryTreatmentVectors.append(temp)
    monthlyBinaryTreatmentVectors = np.array(monthlyBinaryTreatmentVectors)
    binaryTreatmentVectors = []
    # If the aggregate option is chosen, then logically OR each binary string in a segment to form one
    # binary string of length noBits per each segment
    if(aggregateTreatmentsOption):
        boolMonthlyTreatVectors = monthlyBinaryTreatmentVectors.astype(bool)
        for i in range(monthlyBinaryTreatmentVectors.shape[0]):
            rowOfMonths = np.zeros(noBits + signBitOption).astype(bool)
            for j in range(monthlyBinaryTreatmentVectors.shape[1]):
                rowOfMonths = np.logical_or(rowOfMonths, boolMonthlyTreatVectors[i, j])
            binaryTreatmentVectors.append(rowOfMonths.astype(int))
    # If the aggregate treatment option is not chosen, then monthly granujlarity is preserved by
    # Concatenating every binary string in a segment to form one long binary string of length
    # noBits * segmentLength per each segment
    else:
        for i in range(monthlyBinaryTreatmentVectors.shape[0]):
            rowOfMonths = np.array([])
            for j in range(monthlyBinaryTreatmentVectors.shape[1]):
                rowOfMonths = np.concatenate([rowOfMonths, monthlyBinaryTreatmentVectors[i, j]])
            binaryTreatmentVectors.append(rowOfMonths)
    binaryTreatmentVectors = np.array(binaryTreatmentVectors)
    # Initialize distance matrix and compute distances based on the selected distance method
    distanceMatrix = np.zeros((len(binaryTreatmentVectors), len(binaryTreatmentVectors)))
    if binaryDistanceMethod == 1: #Sokal Michener
        for i in range(binaryTreatmentVectors.shape[0]):
            for j in range(i, binaryTreatmentVectors.shape[0]):
                if useManualDistanceMethod:
                    dist = np.sum(np.logical_xor(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))/(noBits + signBitOption)
                else:
                    dist = distance.hamming(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool))
                distanceMatrix[i, j] = dist
                distanceMatrix[j, i] = dist
    elif binaryDistanceMethod == 2: #Jaccard
        for i in range(binaryTreatmentVectors.shape[0]):
            for j in range(i, binaryTreatmentVectors.shape[0]):
                if np.array_equal(binaryTreatmentVectors[i], binaryTreatmentVectors[j]):
                    dist = 0
                else:
                    if useManualDistanceMethod:
                        numerator = np.sum(np.logical_xor(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))
                        denomenator = np.sum(np.logical_or(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))
                        dist = numerator/denomenator
                    else:
                        dist = distance.jaccard(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool))
                distanceMatrix[i, j] = dist
                distanceMatrix[j, i] = dist
    elif binaryDistanceMethod == 3: #Rogers Tanimoto
        for i in range(binaryTreatmentVectors.shape[0]):
            for j in range(i, binaryTreatmentVectors.shape[0]):
                if useManualDistanceMethod:
                    numerator = 2 * np.sum(np.logical_xor(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))
                    denomenator = np.sum(np.logical_xor(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool))) + noBits + signBitOption
                    dist = numerator/denomenator
                else:
                    dist = distance.rogerstanimoto(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool))
                distanceMatrix[i, j] = dist
                distanceMatrix[j, i] = dist
    else: #Sokal Sneath II
        for i in range(binaryTreatmentVectors.shape[0]):
            for j in range(binaryTreatmentVectors.shape[0]):
                if np.array_equal(binaryTreatmentVectors[i], binaryTreatmentVectors[j]):
                    dist = 0
                else:
                    if useManualDistanceMethod:
                        numerator = 2 * np.sum(np.logical_xor(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))
                        denomenator = (2 * np.sum(np.logical_xor(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))) + np.sum(np.logical_and(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool)))
                        dist = numerator/denomenator
                    else:
                        dist = distance.sokalsneath(binaryTreatmentVectors[i].astype(bool), binaryTreatmentVectors[j].astype(bool))
                distanceMatrix[i, j] = dist
    # Normalize distance matrix by dividing the whole thing by the largest value.
    max = np.amax(distanceMatrix)
    distanceMatrix = np.divide(distanceMatrix, max)
    return distanceMatrix

# Converts integer encoded treatments into a binary string of length noBits
def convTreatmentToBinaryArray(treatment, signBitOption, noBits):
    signBit = 0
    if -1 in treatment:
        signBit = 1
    binaryTreatment = np.zeros(noBits)
    if 0 not in treatment:
        listTreatments = list(treatment)
        for i in range(len(listTreatments)):
            if listTreatments[i] != -1:
                binaryTreatment[listTreatments[i] - 1] = 1
    if signBitOption:
        binaryTreatment = np.insert(binaryTreatment, 0, signBit, axis=0)
    return binaryTreatment

# Cluster from a treatment distance matrix using HDB Scan, and print these clusters
# into a txt file
def clusterFromDistMatrix(distanceMatrix):
    clusterer = HDBSCAN(min_cluster_size = 2, metric = 'precomputed')
    clusterer.fit(distanceMatrix)
    labels = clusterer.labels_
    probs = clusterer.probabilities_
    labels = np.array((labels))[np.newaxis]
    labels = labels.T
    probs = np.array((probs))[np.newaxis]
    probs = probs.T
    results = np.concatenate((probs, gv.medSequenceMatrix), axis = 1)
    results = np.concatenate((labels, results), axis = 1)
    results = np.array(sorted(results, key=lambda a_entry: a_entry[0]))
    with open('treatmentClusters.txt', 'w') as csvfile:
        csvfile.write("Cluster; Probability; ID; Month 1; Month 2; Month 3; Month 4; Month 5; Month 6")
        csvfile.write('\n')
        for i in range(results.shape[0]):
            csvfile.write(str(results[i, 0]))
            csvfile.write(';')
            csvfile.write(str(results[i, 1]))
            csvfile.write(';')
            csvfile.write(str(results[i, 2]))
            csvfile.write(';')
            for j in range(3, results.shape[1]):
                csvfile.write(str(results[i, j]).replace('{', '').replace('}', ''))
                csvfile.write(';')
            csvfile.write('\n')
