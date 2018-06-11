#***************************************************************************************#
# Title:        main.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       Executes the various functions that are a part of the capstone
#       project. Allows the user to set the parameters for the execution
#       of said functions. The functions are included in a package called helperFunctions.
#***************************************************************************************#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import datetime
from datetime import datetime, date, time, timedelta
import math
from scipy import stats
from hdbscan import HDBSCAN
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import MinMaxScaler
import time
import inspect
import csv
from ipykernel.tests.test_serialize import point
from sympy.polys.partfrac import apart
from sympy.polys.polytools import intervals
from docutils.writers.docutils_xml import RawXmlError
import os
import glob
import helperFunctions as hf
import globalVariables

#parameters for execution:
alpha = .7 #sets the weighting between lab and treatment values
whichAnalysis = "lab" #decides whether to run lab, treatment, or combinational analysis
whichDist = "P" #decides whether to use pearsons, tau, or spearman's correlation for distance calculations
whichNorm = "M" #chooses a normalization method
overlap = 0 #sets the amount of overlap between bins
segLength = 6 #sets the number of 28-day segments in a bin
graph = 1 #decides whether to use graphing capabilites
noBits = 24 #Set to desired number of bits for the encoded treatment binary string
signBitOption = 1 #Set to 1 if you want to use the sign bit (AKA if you
# want to keep track of unknown drugs by using a sign bit.)
aggregateTreatmentsOption = 0 #Set to 0 for keeping monthly granularity. Set
# to 1 to encode all 6 months as one aggregate binary string
binaryDistanceMethod = 2 #Use to choose the treatment distance calcluation method
# Set binaryDistanceMethod = 1 for Sokal Michener
# Set binaryDistanceMethod = 2 for Jaccard
# Set binaryDistanceMethod = 3 for Rogers Tanimoto
# Set binaryDistanceMethod = 4 for Sokal Sneath II
useManualDistanceMethod = 1
# Set useManualDistanceMethod = 1 to compute distances manually, set useManualDistanceMethod = 0
# to use the SciPy distance calculations for treatments
if(graph == 1):
    path = glob.glob('./graphs/*')
    for file in path:
        os.remove(file)
globalVariables.treatDict = hf.getTreatments()
hf.extractRawInfo()
hf.rawDelete()
hf.rawBinMaker(segLength, overlap)
toBeNormalized = hf.buildFLCMatrix(segLength)
hf.getClustersOnTrend(toBeNormalized, whichDist, whichNorm, segLength, graph)
treatDistMatrix = hf.getDistancesFromMeds(noBits, signBitOption, aggregateTreatmentsOption, binaryDistanceMethod, useManualDistanceMethod, segLength)
hf.clusterFromDistMatrix(treatDistMatrix)
