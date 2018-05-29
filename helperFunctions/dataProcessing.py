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
from datetime import datetime, date, time, timedelta
import globalVariables as gv
import os
import pandas as pd

#Will build a matrix from the raw and interpolated values in
#the specified CSV file. This raw data matrix will then
#be used in our clustering algorithms.
def buildFLCMatrix(segLength):
    with open(gv.labSequenceFileName, "r") as rawPatientData:
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

#Extracts the the treatments from all of the patients in Patient List
#using the LabTreatmentsList csv file. Returns a matrix that will be
#used in the future for building out the treatment dictionary
def getTreatments():
	patList = open("./dataMeasurements/PatientList.csv", "r")
	reader = csv.reader(patList)
	patList = np.array(list(reader))
	patList = np.delete(patList, np.s_[0], axis=0)
	numPats = patList.shape[0]

	treatList = open("./dataMeasurements/LabsTreatmentsList.csv", "r")
	reader = csv.reader(treatList)
	treatList = np.array(list(reader))
	treatList = np.delete(treatList, np.s_[0], axis=0)

	treatmentDictionary = {}
	for i in range(0, treatList.shape[0] - 1):
		if treatList[i][1] != "Lab":
			treatList[i][5] = treatList[i][5].replace(" 00:00:00", "")
			date = datetime.strptime(treatList[i][5], '%Y-%m-%d').date()
			if treatList[i][0] in treatmentDictionary:
				treatmentDictionary[treatList[i][0]].append([treatList[i][2], date])
			else:
				temp = []
				treatmentDictionary[treatList[i][0]] = temp
				treatmentDictionary[treatList[i][0]].append([treatList[i][2], date])

	for i in treatmentDictionary:
		temp = treatmentDictionary[i]
		treatmentDictionary[i] = sorted(temp, key=lambda temp_entry: temp_entry[1])
	return treatmentDictionary

#Builds out a patient dictionary with the raw kappa/lambda values for each
#lab measurement. Each measurement is appropriately timestamped, and these
#will be used to further build the interpolated matrix
def extractRawInfo():
    with open("./dataMeasurements/PatientList.csv", "r") as patList:
        reader = csv.reader(patList)
        patList = np.array(list(reader))
        patDict = {}
        for i in range(1, patList.shape[0]):
            # KEY = MM-# VALUE = whether they are Kappa or Lambda
            patDict[patList[i][0]] = patList[i][3]

    with open("./dataMeasurements/raw_KappaFLC.csv", "r") as rawKappaData:
        reader = csv.reader(rawKappaData)
        rawKapFLC = np.array(list(reader))

    with open("./dataMeasurements/raw_LambdaFLC.csv", "r") as rawLambdaData:
        reader = csv.reader(rawLambdaData)
        rawLamFLC = np.array(list(reader))

    rawKapFLC = np.delete(rawKapFLC, np.s_[0], axis = 0)
    rawLamFLC = np.delete(rawLamFLC, np.s_[0], axis = 0)
    rawKapFLC = np.delete(rawKapFLC, np.s_[3,4], axis = 1)
    rawLamFLC = np.delete(rawLamFLC, np.s_[3,4], axis = 1)

    #Loop through Kappa patient records, deleting any records of Lambda patients
    #At the same time, read the date into a list of datetime objects for better use
    temp1 = 0
    temp2 = 0
    rawKapDates = []
    rawLamDates = []

    while temp1 != rawKapFLC.shape[0]:
        if patDict.get(rawKapFLC[temp1][0]) != "Kappa":
            rawKapFLC = np.delete(rawKapFLC, np.s_[temp1], axis = 0)
        else:
            rawKapFLC[temp1][2] = rawKapFLC[temp1][2].replace(" 00:00:00", "")
            if len(rawKapFLC[temp1][2]) == 8:
                rawKapFLC[temp1][2] = "0" + rawKapFLC[temp1][2]
                rawKapFLC[temp1][2] = rawKapFLC[temp1][2][0:3] + "0" + rawKapFLC[temp1][2][3:]
            elif len(rawKapFLC[temp1][2]) == 9:
                if rawKapFLC[temp1][2][0:2] == "11" or rawKapFLC[temp1][2][0:2] == "12" or rawKapFLC[temp1][2][0:2] == "10":
                    rawKapFLC[temp1][2] = rawKapFLC[temp1][2][0:3] + "0" + rawKapFLC[temp1][2][3:]
                else:
                    rawKapFLC[temp1][2] = "0" + rawKapFLC[temp1][2]
            rawKapDates.append(datetime.strptime(rawKapFLC[temp1][2], '%Y-%m-%d').date())
            temp1 = temp1 + 1

    #Do the same thing for the Lambda patient records
    while temp2 != rawLamFLC.shape[0]:
        if patDict.get(rawLamFLC[temp2][0]) != "Lambda":
            rawLamFLC = np.delete(rawLamFLC, np.s_[temp2], axis = 0)
        else:
            rawLamFLC[temp2][2] = rawLamFLC[temp2][2].replace(" 00:00:00", "")
            rawLamDates.append(datetime.strptime(rawLamFLC[temp2][2], '%Y-%m-%d').date())
            temp2 = temp2 + 1

    rawKapDates = np.array(rawKapDates)
    rawLamDates = np.array(rawLamDates)

    gv.raw_X = np.concatenate((rawKapFLC, rawLamFLC), axis=0)

    gv.raw_dates = np.concatenate((rawKapDates, rawLamDates), axis=0)

    gv.dataTest = {}

    ## PATIENT, FLC VALUE, TEST DATE
    ## N X M MATRIX = ROWS X COLUMNS
    for i in range(0, gv.raw_X.shape[0]):
        # If this patient is already in the dictionary, add additional FLC test values
        if gv.raw_X[i][0] in gv.dataTest:
            gv.dataTest[gv.raw_X[i][0]].append([gv.raw_X[i][1], gv.raw_dates[i]])
        else: # add to dictionary
            gv.dataTest[gv.raw_X[i][0]] = []
            gv.dataTest[gv.raw_X[i][0]].append([gv.raw_X[i][1], gv.raw_dates[i]])
    
    #THESE ARE FOR TESTING. TO BE REMOVED
    with open("rawValuesMatrix.csv", "w") as csvfile:
        for key,value in gv.dataTest.items():
            csvfile.write(key + ': ')
            for v in value:
                csvfile.write(v[0] + ", ")
                csvfile.write(str(v[1]) + " ")
            csvfile.write('\n')

    with open("treatmentDictionary.csv", "w") as csvfile:
        for key,value in gv.treatDict.items():
            csvfile.write(key + ': ')
            for v in value:
                csvfile.write(v[0] + ", ")
                csvfile.write(str(v[1]) + " ")
            csvfile.write('\n')

#deletes patients that had treatments before their baseline test or less
#six months of data. Also builds the smoldering patients dictionary that
#we do not use for now.
def rawDelete():
    for i in gv.dataTest.keys():
        temp = gv.dataTest[i]
        gv.dataTest[i] = sorted(temp, key=lambda temp_entry: temp_entry[1])

    keysToDelete = []
    smolderingRawPatientsDict = {}

    for i in gv.dataTest.keys():
        tempFLC = gv.dataTest[i]
        if((tempFLC[np.array(gv.dataTest[i]).shape[0] - 1][1] - tempFLC[0][1]).days <= 180):
            # print("patient with less than six months: " + i)
            keysToDelete.append(i)
        else:
            if i not in gv.treatDict.keys():
                smolderingRawPatientsDict[i] = gv.dataTest[i]
                keysToDelete.append(i)
                # print("smoldering patient: " + i)
            else:
                tempTreat = gv.treatDict[i]
                if(tempFLC[0][1] > tempTreat[0][1]):  # #tempFLC[row][column] -> FLC DATE KAPPA/LAMBDA*
                    keysToDelete.append(i)
                    # print("patient with treatment before reading: " + i)
                haveFoundSixMonth = False
                for j in range(0, np.array(gv.dataTest[i]).shape[0]):  # for every row in matrix
                    if (((tempFLC[j][1] - tempFLC[0][1]).days >= 180) and (haveFoundSixMonth != True)):
                        sixMonthIndex = j
                        haveFoundSixMonth = True
                firstSixMonths = np.array(tempFLC)[:(sixMonthIndex), :]
                #dataTest[i] = firstSixMonths
                #tempFLC = dataTest[i]
                # print("patient: " + i)
                # else:
                    # print("good patient: " + i)
    for i in keysToDelete:
        del gv.dataTest[i]

    #FOR TESTING. TO BE REMOVED
    with open("rawValuesFilter_1.csv", "w") as csvfile:
        for key,value in gv.dataTest.items():
            csvfile.write(key + ': ')
            for v in value:
                csvfile.write(v[0] + ", ")
                csvfile.write(str(v[1]) + " ")
            csvfile.write('\n')

#Looks at the dictionary data test to create a nested dictionary called outerdict
#the key is the bin number and the values are a dicitonary of items in that list.
#Can input the segment length of each patient and how much we want to overlap them.
def rawBinMaker(segmentLength, overlapBy):
    ## MAP
    ## KEY = MM-
    ## VALUE = DICTIONARY => KEY = NUMBER OF WEEKS - (3->5)*N
    ##                       VALUE = LIST OF FLC VALUES IN THAT TIME FRAME
   
    #Created medSequenceMatrix and labSequenceMatrix which is the csv files in matrix form
    #and without dates. We may not need labSequenceMatrix and we can delete it later if unused.
    gv.medSequenceMatrix = []
    gv.labSequenceMatrix = []
    gv.medSequenceFileName = 'medicationSL_' + str(segmentLength) + '_OL_' + str(overlapBy) + '.csv'
    gv.labSequenceFileName = 'labValuesSL_' + str(segmentLength) + '_OL_' + str(overlapBy) + '.csv'
    try:
        os.remove(gv.medSequenceFileName)
    except OSError:
        pass
    try:
        os.remove(gv.labSequenceFileName)
    except OSError:
        pass
    gv.outerDict = {}
    weekCounter = 1
    for eachKey, value in gv.dataTest.items():
        allTests = []
        allDates = []
        index = 0;
        for v in value:
            if (str(v[1])) in allDates:
                allTests[index - 1] = v[0]
                allDates[index - 1] = str(v[1])
            else:
                allTests.append(v[0])
                allDates.append(str(v[1]))
                index = index + 1
        gv.outerDict[eachKey] = properSampleMaker(eachKey, allTests, allDates, segmentLength, overlapBy)
    with open("rawValuesBins_1.csv", "w", newline='') as csvfile:
        for x in gv.outerDict:
            # X is the Patient ID number
            csvfile.write(x)
            csvfile.write("\n")
            for y in gv.outerDict[x]:
                treatment = str(gv.outerDict[x][y][2]).replace(",", ";")
                temp = str(gv.outerDict[x][y][0]) + ", " +  str(gv.outerDict[x][y][1]) + ", " + treatment
                temp = temp.strip("'[]")
                csvfile.write("Bin: " + str(y) + ", " + temp.replace("'", ""))
                csvfile.write("\n")
    gv.labSequenceMatrix = np.array(gv.labSequenceMatrix)
    gv.medSequenceMatrix = np.array(gv.medSequenceMatrix)

#Makes the inner dictionary for each patient
def properSampleMaker(patientID, FLC_Value, Date, segmentLength, overlapBy):
    dataDict = {'Date': Date, 'FLC_Value': FLC_Value}
    df = pd.DataFrame(dataDict, columns = ['Date', 'FLC_Value'])
    df = df.set_index(pd.DatetimeIndex(df['Date']))
    del df['Date']
    df['FLC_Value'] = df['FLC_Value'].apply(pd.to_numeric, errors='coerce')
    resample = df.resample('D').mean()
    interpolated = resample.interpolate(method='linear')
    downSample = interpolated.resample('28D').first()
    finalData = downSample.reset_index()
    finalData = finalData.values
    innerDict = {}
    datesForTreat = []
    if (patientID == 'MM-120'):
#         with open('MM120_test.csv', 'w', newline='') as csvfile:
#             csvfile.write('\n'.join('{}, {}, {}'.format(resample[x], interpolated[x], downSample[x]) for x in range(0, len(interpolated.values))))
        #print("Original")
        #print(df)
        #print("Daily Upsample")
        #print(resample)
        #print("Interpolating Daily")
        #print(interpolated)
        #print("Down Sample 28 Day Interpolation")
        #print(downSample)
        pass

    for x in range(len(finalData) - 1):
        #treatmentArray = treatmentAdder(finalData[:, 0], np.array(treatDict[patientID]), len(finalData) - 1, patientID)
        innerDict[x] = [finalData[x, 0].to_pydatetime().strftime('%Y-%m-%d'), round(finalData[x, 1], 2)] # Date, FLC
        datesForTreat.append(finalData[x, 0].to_pydatetime().strftime('%Y-%m-%d'))
    gv.treatMatrix = tryTreatment(datesForTreat, patientID)
    indexT = 0
    counterT = 0
    with open(gv.medSequenceFileName, 'a') as csvfile:    #change to file name
        while(indexT + segmentLength < len(datesForTreat)):
            temp = []
            csvfile.write(str(patientID) + "." + str(counterT))
            temp.append(str(patientID) + "." + str(counterT))
            csvfile.write(",")
            for i in range(indexT, segmentLength + indexT):
                csvfile.write(str(gv.treatMatrix[i]).replace(',', ';') + "," + datesForTreat[i])
                temp.append(gv.treatMatrix[i])
                csvfile.write(',')
            csvfile.write('\n')
            gv.medSequenceMatrix.append(temp)
            indexT = indexT + segmentLength - overlapBy
            counterT += 1
    indexL = 0
    counterL = 0
    with open(gv.labSequenceFileName, 'a') as csvfile:  #change to file name
        while(indexL + segmentLength < len(datesForTreat)):
            temp = []
            csvfile.write(str(patientID) + "." + str(counterL))
            temp.append(str(patientID) + "." + str(counterL))
            csvfile.write(",")
            for i in range(indexL, segmentLength + indexL):
                tempList = innerDict[i]
                csvfile.write(str(tempList[1]) + "," + datesForTreat[i])
                temp.append(tempList[1])
                csvfile.write(',')
            csvfile.write('\n')
            gv.labSequenceMatrix.append(temp)
            indexL = indexL + segmentLength - overlapBy
            counterL += 1
    for x in range(len(finalData) - 1):
        innerDict[x].append(gv.treatMatrix[x])
    return innerDict

#builds the treatment array
def tryTreatment(Dates, Patient):
    #print
    # (Dates)
    binDates = []
    treatments = []
    temp = {}
    for i in range(0, len(Dates)):
        binDates.append(datetime.strptime(Dates[i], '%Y-%m-%d').date())
    patientsArray = np.array(gv.treatDict[Patient])
    numData = patientsArray.shape[0]
    temp[0] = set()
    temp[0].add(0)
    for i in range(0, numData):
        treatment = patientsArray[i, 0]
        currDate = patientsArray[i, 1]
        for j in range(1, len(binDates)):
            if currDate >= binDates[j-1] and currDate < binDates[j]:
                if j not in temp.keys():
                    temp[j] = set()
                if(treatment == 'Lenalidomide'):
                    temp[j].add(1)
                elif(treatment == 'Bortezomib'):
                    temp[j].add(2)
                elif(treatment == 'Carfilzomib'):
                    temp[j].add(3)
                elif(treatment == 'Dexamethasone'):
                    temp[j].add(4)
                elif(treatment == 'Pomalidomide'):
                    temp[j].add(5)
                elif(treatment == 'Thalidomide'):
                    temp[j].add(6)
                elif(treatment == 'Cyclophosphamide'):
                    temp[j].add(7)
                elif(treatment == 'Melphalan'):
                    temp[j].add(8)
                elif(treatment == 'Prednisone'):
                    temp[j].add(9)
                elif(treatment == 'Ixazomib'):
                    temp[j].add(10)
                elif(treatment == 'Cisplatin'):
                    temp[j].add(11)
                elif(treatment == 'Doxorubicin'):
                    temp[j].add(12)
                elif(treatment == 'Etoposide'):
                    temp[j].add(13)
                elif(treatment == 'Vincristine'):
                    temp[j].add(14)
                elif(treatment == 'Daratumumab'):
                    temp[j].add(15)
                elif(treatment == 'Elotuzumab'):
                    temp[j].add(16)
                elif(treatment == 'Bendamustine'):
                    temp[j].add(17)
                elif(treatment == 'Panobinostat'):
                    temp[j].add(18)
                elif(treatment == 'Venetoclax'):
                    temp[j].add(19)
                elif(treatment == 'CAR-T'):
                    temp[j].add(20)
                else:
                    temp[j].add(-1)
    for k in range(0, len(binDates)):
        if k not in temp.keys():
            temp[k] = set()
            temp[k].add(0)
    for i in range(0, len(binDates)):
        treatments.append(temp[i])
    return treatments