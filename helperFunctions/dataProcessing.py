#***************************************************************************************#
# Title:        dataProcessing.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Version History:
#       05/15/2018      First Release
#
# Description:
#       Runs the data collection algorithms to prepare for clustering of patients.
#       Requires the following types of files:
#           1. Uninterpolated Kappa and Lambda FLC values for all patients
#           2. List of patients defining type of Myeloma per patient
#           3. CONTINUE
#       Will provide interpolation on the lab values based on overlap and
#       and segment length selection.
#***************************************************************************************#
import csv
import numpy as np

FLCCSVPath = "miniSeqsLabsSL_6_OL_0_BENCHMARK.csv" #path to CSV file containing binned patients

#Will build a matrix from the raw and interpolated values in
#the specified CSV file. This raw data matrix will then
#be used in our clustering algorithms.
def buildFLCMatrix(segLength):
    with open(FLCCSVPath, "r") as rawPatientData:
        reader = csv.reader(rawPatientData)
        FLCMatrix = np.empty(segLength + 1)
        for line in list(reader):
            column = 1
            justData = [] #this matrix will have patient number and each data point
            justData.append(line[0])
            while column < len(line) - 1:
                justData.append(float(line[column].split(",", 1)[0])) #gets rid of the date attached to each FLC value
                column = column + 2
            FLCMatrix = np.vstack((FLCMatrix, justData))
        FLCMatrix = np.delete(FLCMatrix, (0), axis=0)
        return FLCMatrix