#***************************************************************************************#
# Title:        main.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Version History:
#       05/15/2018      First Release
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
whichDist = "E" #decides whether to use pearsons, tau, or spearman's correlation for distance calculations
whichNorm = "N" #chooses a normalization method
overlap = 0 #sets the amount of overlap between bins
segLength = 6 #sets the number of 28-day segments in a bin
graph = 0 #decides whether to use graphing capabilites
globalVariables.treatDict = hf.getTreatments()
hf.extractRawInfo()
hf.rawDelete()
hf.rawBinMaker(segLength, overlap)
toBeNormalized = hf.buildFLCMatrix(segLength)
hf.getClustersOnTrend(toBeNormalized, whichDist, whichNorm)