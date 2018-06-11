#***************************************************************************************#
# Title:        combinationClustering.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       All functions necessary for factoring in both the FLC and treatment values
#       in the clustering. Can change the weighting of treatment:FLC distances by
#       modifying alpha.
#***************************************************************************************#

### pass in treatment and lab distance matrix into this function, returns the combined distance matrix

# matrix1 = TREATMENTS matrix2 = LABS
# overallproj.py parameter input: matrix 1 = medicationDistNorm matrix 2 = pearsonDistNorm
# overallproj.py output = combDistMatrix
def combineDistances(matrix1, matrix2): 
    overallDist = np.zeros((matrix1.shape[0], matrix1.shape[0]))
    if (matrix1.shape[0] == matrix2.shape[0]):
        for i in range(0, matrix1.shape[0]):
            for j in range(0, matrix1.shape[1]):
                overallDist[i][j] = (treatmentAlpha * matrix1[i][j]) + ((1 - treatmentAlpha) * matrix2[i][j])
    np.savetxt("combinedDistanceMatrix.csv", overallDist, delimiter=",")
    return overallDist

### then normalize this combined distance matrix 
# overallproj.py input = combDistMatrix
# overallproj.py output = combDistNorm
def normalizeDist(matrix): 
    oneDimMatrix = np.ndarray.flatten(matrix)
    #print(oneDimMatrix)
    maxValue = np.amax(oneDimMatrix)
    print("Maximum from matrix: " + str(maxValue))
    maxValue = max(maxValue, denominatorTolerance)
    print("Final max value: " + str(maxValue))
#     print(matrix)
#     print()
#     print(max)
#     print()
    newMatrix = np.divide(matrix, maxValue)
    #print(newMatrix)
    return newMatrix
  
  
### then perform a clustering on the combined normalize distance matrix 
# overallproj.py input: distanceMatrix = combDistNorm, selector = "combinedNormalizedDistanceClusters"
def clusterFromDistMatrix(distanceMatrix, selector):
    clusterer = HDBSCAN(min_cluster_size = 2, metric = 'precomputed')
    clusterer.fit(distanceMatrix)
    labels = clusterer.labels_
    probs = clusterer.probabilities_
    labels = np.array((labels))[np.newaxis]
    labels = labels.T
    probs = np.array((probs))[np.newaxis]
    probs = probs.T
    if (selector == "treatmentClusters"):
        results = np.concatenate((probs, medSequenceMatrix[0:50, :]), axis = 1)
        results = np.concatenate((labels, results), axis = 1)
        results = np.array(sorted(results, key=lambda a_entry: a_entry[0]))
    elif (selector == "combinedNormalizedDistanceClusters"):
        #print(benchmarkLabs[:, 0].T)
        results = np.column_stack((probs, benchmarkLabs[:, 0]))
        results = np.concatenate((labels, results), axis = 1) 
        
    results = np.array(sorted(results, key=lambda a_entry: a_entry[0]))   
        
    if (selector == "combinedNormalizedDistanceClusters"): 
        #print(results)
        # hdb ___ plotting ###
        #hdb_unique_labels = set(hdb_labels)
        hdb_colors = plt.cm.Spectral(np.linspace(0, 1, len(labels)))
        fig = plt.figure(figsize=plt.figaspect(0.5))
        hdb_axis = fig.add_subplot('121')
     
        for i, (k, col) in enumerate(zip(labels, hdb_colors)):
            if k == -1:
                # Black used for noise.
                col = 'k'
             # Not sure what the x, y axes should be, this is how I had it set up for the previous run of HDBSCAN
             # Right now it's only comparing the first and second data points of all patients
             
            print(i)
            print(results[i, :])
            print(benchmarkLabs[i, 1:])
#             hdb_axis.plot([1, 2, 3, 4, 5, 6, 7], clusterer[i, :].astype(float), 'o', markerfacecolor=col, markeredgecolor='k', markersize=6)
#         hdb_axis.set_title('HDBSCAN\nEstimated number of clusters: %d' % n_clusters_hdb_)
#         plt.show()
        
#         for clusterRow in range(0, clusterMatrix.shape[0]): #plotting
#                  plt.plot([1,2,3,4,5,6], clusterMatrix[clusterRow,1:].astype(float), label=clusterMatrix[clusterRow, 0])
#                  plt.legend()
        
        with open(selector + '.csv', 'w') as csvfile:
            csvfile.write("Cluster; Probability; ID; Month 1; Month 2; Month 3; Month 4; Month 5; Month 6")
            csvfile.write('\n')
            for i in range(results.shape[0]):
                csvfile.write(str(results[i, 0]))
                csvfile.write(',')
                csvfile.write(str(results[i, 1]))
                csvfile.write(',')
                csvfile.write(str(results[i, 2]))
                csvfile.write('\n')
    if (selector == "treatmentClusters"):
        with open(selector + '.csv', 'w') as csvfile:
            csvfile.write("Cluster; Probability; ID; Month 1; Month 2; Month 3; Month 4; Month 5; Month 6")
            csvfile.write('\n')
            for i in range(results.shape[0]):
                csvfile.write(str(results[i, 0]))
                csvfile.write(',')
                csvfile.write(str(results[i, 1]))
                csvfile.write(',')
                csvfile.write(str(results[i, 2]))
                csvfile.write(',')
                for j in range(3, results.shape[1]):
                    csvfile.write(str(results[i, j]).replace('{', '').replace('}', ''))
                    csvfile.write(',')
                csvfile.write('\n')
