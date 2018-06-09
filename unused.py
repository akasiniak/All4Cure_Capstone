#***************************************************************************************#
# Title:        unused.py
# Authors:      Alexander Kasiniak
#               Cece Landau
#               Kevin Lau
#
# Description:
#       List of unused functions that could be of use at some point in the
#       future.
#***************************************************************************************#

def derivativeMaker():
    # Create first derivative column
    D1 = np.zeros((len(X), 1))
    for i in range(1, len(dates) - 1):
        if X[i][0] == X[i - 1][0]:
            xdif = dates[i] - dates[i - 1]
            ydif = np.float(X[i][1]) - np.float(X[i - 1][1])
            if ydif == 0:
                D1[i] = 0
            elif xdif.total_seconds() / 86400 == 0:
                if ydif > 0:
                    D1[i] = float('inf')
                else:
                    D1[i] = float('-inf')
            else:
                D1[i] = ydif / (xdif.total_seconds() / 86400)
    min_val = D1[np.isfinite(D1)].min()
    max_val = D1[np.isfinite(D1)].max()
    for i in range(1, len(dates) - 1):
        if D1[i] == float('-inf'):
            D1[i] = min_val
        elif D1[i] == float('inf'):
            D1[i] = max_val

    # Create second derivative column
    D2 = np.zeros((len(X), 1))
    for i in range(1, len(dates) - 1):
        if X[i][0] == X[i - 1][0]:
            xdif = dates[i] - dates[i - 1]
            ydif = np.float(D1[i]) - np.float(D1[i - 1])
            if ydif == 0:
                D2[i] = 0
            elif xdif.total_seconds() / 86400 == 0:
                if ydif > 0:
                    D2[i] = float('inf')
                else:
                    D2[i] = float('-inf')
            else:
                D2[i] = ydif / (xdif.total_seconds() / 86400)
    min_val2 = D2[np.isfinite(D2)].min()
    max_val2 = D2[np.isfinite(D2)].max()
    for i in range(1, len(dates) - 1):
        if D2[i] == float('-inf'):
            D2[i] = min_val2
        elif D2[i] == float('inf'):
            D2[i] = max_val2

    return D1, D2

def patientBinCreator(patientID, FLC_Value, Date):
    #innerDict = {'0' : ['test initial'], '1' : ['a', 'b', 'c']}
    for x in range(0, len(Date)):
        Date[x] = datetime.strptime(Date[x], '%Y-%m-%d').date()
    innerDict = {}
    innerDict[0] = [FLC_Value[0], Date[0]]
    binNumber = 0
    dateToCompare = Date[0]
    for dateInstance in range(1, len(Date) - 1):
        timeDiff = Date[dateInstance] - dateToCompare
        timeDifference = timeDiff.days
        # If it's within 21 to 35 days from the previous date marker, put in current bin
        if (timeDifference > 21 and timeDifference < 35):
            if (binNumber == 0):
                binNumber = 1
            if (binNumber in innerDict):
                innerDict[binNumber].extend([FLC_Value[dateInstance], Date[dateInstance]])
            else:
                innerDict[binNumber] = [FLC_Value[dateInstance], Date[dateInstance]]
        # Time difference between current point and date to compare falls outside of the bin range,
        # increase the bin number by one, and later fill that bin with an interpolated (linear) value
        elif (timeDifference > 35):
            binNumber = binNumber + 1
            #print("Not within range, starting bin number " + str(binNumber) + "!")
            treatmentsArray = treatDict[patientID]
            currentTreat = []
            oldTreat = []
            for x in range(0, len(treatmentsArray) - 1):
                # Wait to actually update dateToCompare so we can easily fill oldTreat array
                dateInterest = treatmentsArray[x][1] - dateToCompare + timedelta(days=28)
                dateInt = dateInterest.days
                if (dateInt > 21 and dateInt < 35):
                    if (treatmentsArray[x][0] not in currentTreat):
                        currentTreat.append(treatmentsArray[x][0])
                if (timeDifference > 21 and timeDifference < 35):
                    if (treatmentsArray[x][0] not in oldTreat):
                        oldTreat.append(treatmentsArray[x][0])
            dateToCompare = dateToCompare + timedelta(days=28)
            if (binNumber in innerDict):
                innerDict[binNumber].extend(["-1", dateToCompare])
            else:
                innerDict[binNumber] = ["-1", dateToCompare]
        else:
#             print("Too small for range, delete this item")
#             print(str(timeDifference))
            pass
    print(innerDict)
    print("About to return to rawBinMaker")
    return innerDict

def FLCdictionary(D1, D2):
    global FLCdict
    FLCdict = {}

    for i in range(0, X.shape[0]):
    	if X[i][0] in FLCdict:
    		FLCdict[X[i][0]].append([X[i][1], dates[i], D1[i], D2[i]])
    	else:
    		FLCdict[X[i][0]] = []
    		FLCdict[X[i][0]].append([X[i][1], dates[i], D1[i], D2[i]])


    smolderingPatientsDict = {}

    for i in FLCdict.keys():
        temp = FLCdict[i]
        FLCdict[i] = sorted(temp, key=lambda temp_entry: temp_entry[1])

    keysToDelete = []
    for i in FLCdict.keys():
        tempFLC = FLCdict[i]
        if((tempFLC[np.array(FLCdict[i]).shape[0] - 1][1] - tempFLC[0][1]).days <= 180):
            # print("patient with less than six months: " + i)
            keysToDelete.append(i)
        else:
            if i not in treatDict.keys():
                smolderingPatientsDict[i] = FLCdict[i]
                keysToDelete.append(i)
                # print("smoldering patient: " + i)
            else:
                tempTreat = treatDict[i]
                if(tempFLC[0][1] > tempTreat[0][1]):  # #tempFLC[row][column] -> FLC DATE KAPPA/LAMBDA*
                    keysToDelete.append(i)
                    # print("patient with treatment before reading: " + i)
                haveFoundSixMonth = False
                for j in range(0, np.array(FLCdict[i]).shape[0]):  # for every row in matrix
                    if (((tempFLC[j][1] - tempFLC[0][1]).days >= 180) and (haveFoundSixMonth != True)):
                        sixMonthIndex = j
                        haveFoundSixMonth = True
                firstSixMonths = np.array(tempFLC)[:(sixMonthIndex), :]
                FLCdict[i] = firstSixMonths
                tempFLC = FLCdict[i]
                # print("patient: " + i)
                # else:
                    # print("good patient: " + i)
    for i in keysToDelete:
        del FLCdict[i]
    ### plotting flc value for each patient ###
    # for i in FLCdict.keys():
    #    tempFLC = np.array(FLCdict[i])
    #    plt.figure()
    #    plt.plot(tempFLC[:, 1], tempFLC[:, 0])
    #    plt.title(i)
    #    # print(tempFLC[:, 0])
    #    # print(tempFLC[:, 1])
    #    plt.show()
    return

def processingWrite():
    global preSpearman
    preSpearman = []
    global useablePatients
    useablePatients = []
    lengthSegment = 5  # # Should we leave this as different than numberOfPoints?
    with open('processed.csv', 'w') as csvfile:
        temp = []
        for i in FLCdict.keys():
            temp = np.array(FLCdict[i])
            numReadings = temp.shape[0]  # ##numReadings = number of readings left to process
            counter = 0
            while(numReadings >= lengthSegment):  # will need to change 5 to a field
                temp2 = []
                for j in range((lengthSegment - 1) * counter, lengthSegment + (lengthSegment - 1) * counter):
                    temp2.append(temp[j][0])
                csvfile.write(str(i) + "-" + str(counter))
                useablePatients.append(str(i) + "-" + str(counter))
                for j in range (((lengthSegment - 1) * counter), lengthSegment + ((lengthSegment - 1) * counter)):
                    csvfile.write(", " + str(temp[j][0]))
                    csvfile.write(", " + str(temp[j][1]))
                csvfile.write('\n')
                preSpearman.append(temp2)
                numReadings = numReadings - lengthSegment + 1
                counter += 1
    preSpearman = np.array(preSpearman)
    preSpearmanNum = np.array(preSpearman.astype(float))

    unprocessedMatrix = np.array(list(zip(useablePatients, np.array(preSpearman.astype(float)))), dtype=object)
    with open('unscaledData.csv', 'w') as csvfile:
        csvfile.write("Patient Number + UnScaled FLC Values" + '\n')
        csvfile.write('\n'.join('{}, {}'.format(x[0], x[1]) for x in unprocessedMatrix))
    return

def extractInfo():
    unpack = getVectors()
    kapFLC = unpack[0]
    kapDates = unpack[2]
    lamFLC = unpack[3]
    lamDates = unpack[5]
    global X
    X = np.concatenate((kapFLC, lamFLC), axis=0)
    global dates
    dates = np.concatenate((kapDates, lamDates), axis=0)
    global treatDict
    treatDict = getTreatments()
    return