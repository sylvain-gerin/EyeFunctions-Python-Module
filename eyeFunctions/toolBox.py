
""" this module provides a series of functions used to handle and operate on python lists of various dimensions
Dependencies: matplotlib.pyplot, scipy.signal, math, csv
Author: Sylvain Gerin, sylvain.gerin@uclouvain.be 
"""

import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from math import atan2, degrees
import csv

def unpack(list2D):
    """ change a 2D list into a 1D one"""
    unpackedList = []
    for thisList in list2D:
        unpackedList += thisList
    return unpackedList

def loadFile(fileName, delimiter='\t', encoding=None):
    """ open a file with csv.reader using a chosen delimiter
    arguments:
    fileName -- the name of the file to open
    delimiter -- the column separator. Default '\t'
    encoding -- the name of the encoding used to decode the file
    """
    with open(fileName, newline ='', encoding=encoding) as file:
        data = csv.reader(file, delimiter=delimiter)
        data = list(data)
    return data

def mean(list1D):
    """ return the mean of a 1D list of float or integers. Omit Nan"""
    list1D = [i for i in list1D if str(i) != 'nan']
    total = sum(list1D)
    try:
        meanValue = total/len(list1D)
    except:
        meanValue = float('nan')
    
    return meanValue

def median(list1D):
    """ return the median of a 1D list of float or integers"""
    list1D = [i for i in list1D if str(i) != 'nan']
    list1D.sort()
    
    n = len(list1D)
    if n % 2 != 0:
        indexToLook = int((n+1)/2) - 1 # python starts at 0
        medianValue = list1D[indexToLook]
    else:
        indexToLook = [int(n/2) - 1, int((n/2)+1) - 1]
        valuesToLook = list1D[indexToLook[0]:indexToLook[1] + 1]
        medianValue = mean(valuesToLook)
    
    return medianValue


def variance(list1D):
    """ return the variance of a list of float or integers"""
    list1D = [i for i in list1D if str(i) != 'nan']
    meanValue = mean(list1D)
    squaredDiff = 0.0
    for i in list1D:
        squaredDiff += (i - meanValue)**2
    variance = squaredDiff / (len(list1D)-1)
    return(variance)

def getMAD(list1D, coefficient=1.4826):
    """ return the median absolute deviation of a list of values and multiply it by a chosen coefficient (default 1.4826)"""
    list1D = [i for i in list1D if str(i) != 'nan']
    toSubtract = median(list1D)
    absoluteDeviations = [abs(i - toSubtract) for i in list1D]
    MAD = median(absoluteDeviations)
    MAD *= coefficient
    return MAD


def sqrt(x):
    """ return the square root of a given number."""
    return(x**0.5)

def nbOfDimensions(dataSet):
    """ inform if a list has 0, 1, 2 or 3 dimensions"""
    if type(dataSet) is not list:
        nbOfLevels = 0
    
    elif type(dataSet[0]) is list:
        if type(dataSet[0][0]) is list:
            nbOfLevels = 3
        else:
            nbOfLevels = 2
    else:
        nbOfLevels = 1
    return nbOfLevels

def stringToFloat(dataSet, zeroIsNan=False):
    """ convert a 1D or 2D list of values into float. If not possible, replace the value with nan
    arguments:
    dataSet -- a 1D or 2D list of strings
    zeroIsNan -- decide whether '0.0' and '0' values are changed into Nan or not. Default False
    """  
    nbOfLevels = nbOfDimensions(dataSet)
    floatedDataset = []
    
    if nbOfLevels == 2:
        for thisTrial in dataSet:
            thisTrialFloated = []
            for thisFrame in range(len(thisTrial)):
                try:
                    if float(thisTrial[thisFrame]) == float(0) and zeroIsNan==True:
                        thisTrialFloated[thisFrame] = float('nan')
                    else:
                        thisTrialFloated += [float(thisTrial[thisFrame])]                   
                except:
                    thisTrialFloated += [float('nan')]
            floatedDataset += [thisTrialFloated]
    
    elif nbOfLevels == 1:
        for thisFrame in range(len(dataSet)):
            try:
                if float(dataSet[thisFrame]) == float(0) and zeroIsNan == True:
                    floatedDataset[thisFrame] = float('nan')
                else:
                    floatedDataset += [float(dataSet[thisFrame])]    
            except:
                floatedDataset += [float('nan')]
                
    return floatedDataset

def stringToInt(dataSet, zeroIsNan=False):
    """ convert a 1D or 2D list of values into int. If not possible, replace the value with nan
    arguments:
    dataSet -- a 1D or 2D list of strings
    zeroIsNan -- decide whether '0.0' and '0' values are changed into Nan or not. Default False
    """  
    nbOfLevels = nbOfDimensions(dataSet)
    intDataset = []
    
    if nbOfLevels == 2:
        for thisTrial in dataSet:
            thisTrialInted = []
            for thisFrame in range(len(thisTrial)):
                try:
                    if int(thisTrial[thisFrame]) == int(0) and zeroIsNan==True:
                        thisTrialInted[thisFrame] = float('nan')
                    else:
                        thisTrialInted += [int(thisTrial[thisFrame])]                   
                except:
                    thisTrialInted += [float('nan')]
            intDataset += [thisTrialInted]
    
    elif nbOfLevels == 1:
        for thisFrame in range(len(dataSet)):
            try:
                if int(dataSet[thisFrame]) == int(0) and zeroIsNan == True:
                    intDataset[thisFrame] = float('nan')
                else:
                    intDataset += [int(dataSet[thisFrame])]    
            except:
                intDataset += [float('nan')]
                
    return intDataset

def isOutlier(list1D, maxSD, messageIn='notOutlier', messageOut='outlier'):
    """ determines whether elements of a list are outliers or not, as defined by a chosen number of standard deviations from the mean
    arguments:
    list1D -- the 1D list of all values to consider
    maxSD -- the maximal number of standrad deviations before considering values as outliers
    messageIn -- the message returned by the function if an item is not an outlier
    messageOut -- the message returned by the function if an item is an outlier
    """
    meanBaseline = mean(list1D)
    sdBaseline = sqrt(variance(list1D))
    upLimit = meanBaseline + maxSD*sdBaseline
    downLimit = meanBaseline - maxSD*sdBaseline
    outlierTrial = []
    
    for thisTrial in range(len(list1D)):
        if list1D[thisTrial] > upLimit or list1D[thisTrial] < downLimit:
            outlierTrial += [messageOut]
        else:
            outlierTrial += [messageIn]
    
    return(outlierTrial)

def getColumnIndex(list2D, headerName):
    try:
        index = list2D[0].index(headerName)
        return index
    except:
        print(f'{headerName} not in list')

def saveList(dataSet, fileName, charSep='\t', blockSep='\n', mode='w'):
    """ Save elements of a list of 1 or 2 dimensions in a specified file.
    arguments:
    dataSet -- the 1 or 2D list to save
    fileName -- the name given to the file in which the list is written
    charSep -- the separator between the list elements. Default '\t'
    blockSep -- the separator between several lists written on the same file. Default '\n'
    mode -- the mode of writing in the created file. 'a' to append new elements to the existing file, 'w' to overwrite. Default 'w'
    """
    # Open a file in 'append' or 'write' mode
    output = open(fileName, mode)
    nbOfLevels = nbOfDimensions(dataSet)
    # Loop through each element of the list
    if nbOfLevels == 2:
        for i in range(len(dataSet)):
            for j in range(len(dataSet[i])):
                # Add these elements in the file
                if j < (len(dataSet[i]) - 1):
                    output.writelines(str(dataSet[i][j]) + charSep)
                elif j >= (len(dataSet[i]) - 1):
                    output.writelines(str(dataSet[i][j]))
            output.writelines(blockSep)
            
    elif nbOfLevels == 1:
        for i in range(len(dataSet)):
            if i < (len(dataSet) - 1):
                output.writelines(str(dataSet[i]) + charSep)
            elif i >= (len(dataSet)) - 1:
                output.writelines(str(dataSet[i]))
            output.writelines('\n')
        output.writelines(blockSep)
    else:
        print('invalid number of Dimensions of the input list')
    
    # Close file
    output.close()
    return None

def getValues(dataSet, variableOfInterest):
    """ extracts all values that have the same index in a 2D or 3D list. Could be seen as 'extracting a column'
    arguments:
    dataSet -- the 2D or 3D list in which values are taken. If the list has 2 Dimensions, it will extract the value of the given index in each list
            If the list has 3 dimensions, it will return a 2D list, each sublist being containing the given index of each sublist
    variableOfInterest -- the index for which values will be extracted
    """
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 3:
        getSamples = []
        for thisTrial in dataSet:
            thisTrialSamples = [thisFrame[variableOfInterest] for thisFrame in thisTrial]
            getSamples += [thisTrialSamples]
    
    elif nbOfLevels == 2:
        getSamples = [thisFrame[variableOfInterest] for thisFrame in dataSet]
        
    return(getSamples)

def selectPeriod(dataSet, startPeriod, endPeriod):
    """ select the items of a 1D or 2D list in a specified range
    arguments:
    dataSet -- a 1D or 2D list in which the wanted items will be extracted
    startPeriod -- the index ofthe first item of interest
    endPeriod -- the index of the last item of interest
    
    if the given dataSet is a 1D list, this is similar to a basic list indexing (myList[start:stop])
    If the list is a 2D list, the values of each separate list will be extracted, and it will return a 2D list
    """
    interestStartIndex = startPeriod
    interestEndIndex = endPeriod
    
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 2:
        # copy the samples of interest
        epuredTrials = []
        for thisTrial in dataSet:
            epuredTrials += [thisTrial[interestStartIndex:interestEndIndex]]
    
    elif nbOfLevels == 1:
        epuredTrials = dataSet[interestStartIndex:interestEndIndex]
    
    return epuredTrials

def derivative(list1D):
    """ compute the difference between the ith+1 and the ith value of a list. Return a list of the same length of the initial list, starting with 0"""
    der = [0] + [(list1D[i+1]-list1D[i]) for i in range(len(list1D)-1)]
    return der


def velocitySearch(list1D, rangeOfSearch, eventVelocityThreshold, eventContinuousFrames, stopAtFirst=True, around=False):
    """ iterate in a given range of a list of floats or integers to find an event based on the difference of consecutive values, and return the index at which an event occured.
    If no event has been found, return None
    arguments:
    list1D -- a 1D list to iterate in
    rangeOfSearch -- range(onset, offset) in which an event is looked for
    eventVelocityThreshold -- the minimal difference of values between two consecutive iterations. Negative threshold indicate a decrease
    eventContinuousFrames -- the number of continuous iterations above threshold necessary to consider an event. 1=2 consecutive iterations
    stopAtFirst -- if True, returns the first iteration at which the threshold is reached, if False, returns the last iteration found, default True
    around -- if velocity must stay within a range of +- the given threshold. If True, make sure the threshold is a positive value. default False
    """
    counter = 0
    thisEvent = None
    
    # loop through all items of a list
    for thisFrame in rangeOfSearch:
        # check if the difference between item n and n+1 is superior to a given threshold
        
        if around == True: # if the difference of consecutive values must stay within a limit
            if hasattr(eventVelocityThreshold, '__iter__'):
                isBetween = min(eventVelocityThreshold) <= list1D[thisFrame+1] - list1D[thisFrame] <= max(eventVelocityThreshold)
            else:
                isBetween = min(eventVelocityThreshold, -eventVelocityThreshold) <= list1D[thisFrame+1] - list1D[thisFrame] <= max(eventVelocityThreshold, -eventVelocityThreshold)
                
            if isBetween:
                counter += 1
                # check if the difference has been above threshold for the wanted number of consecutive frames
                if counter >= eventContinuousFrames:
                    if stopAtFirst == True:
                        thisEvent = thisFrame - (counter-1)
                        break
                    elif stopAtFirst == False:
                        if not isBetween:
                            thisEvent = thisFrame +1
                            break
                        else:
                            thisEvent = thisFrame +1
            # if the velocity is below threshold, counter is reset
            else:
                counter = 0
        
        # if the threshold is negative, check for values below the threshold
        elif eventVelocityThreshold < 0 and around == False:
            if list1D[thisFrame+1] - list1D[thisFrame] <= eventVelocityThreshold:
                counter += 1
                # check if velocity has been above threshold for the wanted number of consecutive frames
                if counter >= eventContinuousFrames:
                    if stopAtFirst == True:
                        thisEvent = thisFrame - (counter-1)
                        break
                    elif stopAtFirst == False:
                        # if the threshold is reached and will be lost at the next iteration, take this event and stop search
                        if list1D[thisFrame+2] - list1D[thisFrame+1] > eventVelocityThreshold:
                            thisEvent = thisFrame +1
                            break
                        else:
                            thisEvent = thisFrame +1
            # if the velocity is below threshold, counter is reset
            else:
                counter = 0
        
        elif eventVelocityThreshold > 0 and around == False:
            if around == False:
                if list1D[thisFrame+1] - list1D[thisFrame] >= eventVelocityThreshold:
                    counter += 1
                    # check if velocity has been above threshold for the wanted number of consecutive frames
                    if counter >= eventContinuousFrames:
                        if stopAtFirst == True:
                            thisEvent = thisFrame - (counter-1)
                            break
                        elif stopAtFirst == False:
                        # if the threshold is reached and will be lost at the next iteration, take this event and stop search
                            if list1D[thisFrame+2] - list1D[thisFrame+1] < eventVelocityThreshold:
                                thisEvent = thisFrame +1
                                break
                            else:
                                thisEvent = thisFrame +1
                # if the velocity is below threshold, counter is reset
                else:
                    counter = 0
    
    return thisEvent

def createDataBase(dataSet, nbOfTrials):
    """ combine the information of different 1D and 2D lists in a same list
    arguments:
    dataSet -- a list of 1D or 2D lists, single elements, etc
    nbOfTrials -- the desired final number of rows """
    
    # create a list that will accumulate the information of each given list
    tempDataSet = []
    
    for thisList in dataSet:
        # if not a list but a single value, repeat it to create a list of it
        if type(thisList) is not list:
            list1D = []
            for i in range(nbOfTrials):
                list1D += [thisList]
            tempDataSet += [list1D]
        else:
            # take the entire 1D or 2D list
            tempDataSet += [thisList]
    
    # for each sublist, take the values specific to 1 trial together 
    thisTrial = 0
    generalDataSet = []
    while thisTrial < nbOfTrials:
        thisTrialData = []
        for thisVariable in tempDataSet:
            if type(thisVariable[thisTrial]) is not list:
                thisTrialData += [thisVariable[thisTrial]]
            else:
                thisTrialData += [i for i in thisVariable[thisTrial]]
        generalDataSet += [thisTrialData]
        thisTrial += 1
    
    return generalDataSet

def descriptives(dataSet):
    """ return the mean, SD and SE of a 1D or 2D list"""
    
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 1:
        dataSetMean = mean(dataSet)
        dataSetSD = sqrt(variance(dataSet))
        dataSetSE = dataSetSD/sqrt(len(dataSet))
    
    elif nbOfLevels == 2:
        dataSetMean = []
        dataSetSD = []
        dataSetSE = []
        
        for thisList in dataSet:
            thisListMean = mean(thisList)
            thisListSD = sqrt(variance(thisList))
            thisListSE = sqrt(variance(thisList)) / sqrt(len(thisList))
            
            dataSetMean += [thisListMean]
            dataSetSD += [thisListSD]
            dataSetSE += [thisListSE]
    
    return(dataSetMean, dataSetSD, dataSetSE)

def operateLists(list1, list2, kind):
    """ perform basic arithmetic operations on lists of the same length
    arguments:
    list1 -- the list containing the first operands
    list2 -- the list containing the second operands
    kind -- the type of operation to perform: 'subtract', 'add', 'multiply', 'divide'
    """
    if len(list1) == len(list2):
        if kind == 'subtract':
            finalList = [list1[i] - list2[i] for i in range(len(list1))]
        elif kind == 'add':
            finalList = [list1[i] + list2[i] for i in range(len(list1))]
        elif kind == 'multiply':
            finalList = [list1[i] * list2[i] for i in range(len(list1))]
        elif kind == 'divide':
            finalList = [list1[i] / list2[i] for i in range(len(list1))]
    else:
        print('lists must have the same length')
    return finalList

def cluster(list2D, col, cond):
    """ Group all the rows of a 2D list for which a given column has the same value and return a list of rows (2D list).
    arguments:
    list2D -- the 2D list to search in
    col -- the column index in which the condition applies
    cond -- the clustering criterion
    """
    # Create an empty list
    clusteredList = []
    # Fill it with all the rows responding to the condition
    for i in range(len(list2D)):
        if str(list2D[i][col]) == str(cond):
            clusteredList += [list2D[i]]
        else:
            continue
    # Return the filled list
    return clusteredList

def centerList(dataSet, baselineValues=None, firstAvailable=True):
    """center lists by subtracting specified values from each value of a 1D or 2D list
    arguments:
    dataSet -- the 1D or 2D list to be centered
    baselineValues -- the values to subtract from lists. None, int, float or 1D list. If None, each list will be centered on its firts value. default None
    firstAvailable -- if centering is made from the first value, ensures that it is not a nan. default True
    """
    
    centerFirst = True if baselineValues is None else False
    
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 2:
        if centerFirst == True:
            if firstAvailable == False:
                baselineValues = [i[0] for i in dataSet]
            else:
                baselineValues = []
                for thisTrial in dataSet:
                    if all(str(values) == 'nan' for values in thisTrial):
                        baselineValues += [float('nan')]
                    else:
                        for thisValue in thisTrial:
                            if str(thisValue) != 'nan':
                                baselineValues += [thisValue]
                                break
        else:
            if nbOfDimensions(baselineValues) == 0:
                baselineValues = [baselineValues for i in range(len(dataSet))]
        
        centeredList = []
        for thisBaseline, thisList in zip(baselineValues, dataSet):
            centeredList += [[i - thisBaseline for i in thisList]]
    
    elif nbOfLevels == 1:
        if centerFirst == True:
            if firstAvailable == False:
                baselineValues = dataSet[0]
            else:
                if all(str(values) == 'nan' for values in dataSet):
                    baselineValues = float('nan')
                else:
                    for thisValue in dataSet:
                        if str(thisValue) != 'nan':
                                baselineValues = thisValue
                                break
        
        centeredList = [i - baselineValues for i in dataSet]
    
    return centeredList

def centerFirstFrame(dataSet):
    """subtract values to be expressed with respect to the first item of a list. Works with 1D and 2D lists"""
    nbOfLevels = nbOfDimensions(dataSet)
    
    firstFrameIndex = 0
    if nbOfLevels == 2:
        for thisTrial in range(len(dataSet)):
            toSubtract = dataSet[thisTrial][firstFrameIndex]
            for thisFrame in range(len(dataSet[thisTrial])):
                dataSet[thisTrial][thisFrame] -= toSubtract
    
    elif nbOfLevels == 1:
        toSubtract = dataSet[firstFrameIndex]
        for thisFrame in range(len(dataSet)):
            dataSet[thisFrame] -= toSubtract
    
    return None

def centerBaseline(dataSet, baselineValues):
    """ subtract a previously computed value from a dataset
    arguments:
    dataSet -- the dataset containing the values to center (1D or 2D list)
    baselineValues -- the baseline values to remove from the dataSet
    !!! the variable containing the baseline values must have a dimension less than the dataSet!!!
    """
    nbOfLevels = nbOfDimensions(dataSet)
    
    centeredTrials = []
    
    if nbOfLevels == 2:
        for thisTrial in range(len(baselineValues)):
            toSubtract = baselineValues[thisTrial]
            thisCenteredTrial = []
            for thisFrame in dataSet[thisTrial]:
                thisCenteredTrial += [thisFrame - toSubtract]
            centeredTrials += [thisCenteredTrial]
            
    elif nbOfLevels == 1:
        toSubtract = baselineValues
        for thisFrame in dataSet:
            centeredTrials += [thisFrame-toSubtract]
    
    return centeredTrials

def plotTrials(list2D, filename, legend=None, fill=None, color=None, show=False):
    """ plot the values of 2D lists, 1 line per list
    arguments:
    list2D -- the 2D list containing each value to show on the same plot
    fileName -- the name under which the plot must be saved
    legend -- a 1D list of labels corresponding to every sublist of the list2D. Default None
    fill -- a 2D list of pairs corresponding to the upper and lower limits of shaded area (for example, mean SE). Default None
    color -- a 1D list of colors, one per sublist of the list2D. Default None
    show -- display the plot. default False
    """
    a = []
    for thisTrial in range(len(list2D)):
        if color is not None:
            a += [plt.plot(list2D[thisTrial], color=color[thisTrial])]
        else:
            a += [plt.plot(list2D[thisTrial], color=color)]
        
        if fill is not None and color is not None:
            plt.fill_between(range(len(list2D[thisTrial])), fill[thisTrial][0], fill[thisTrial][1], alpha = 0.2, color=color[thisTrial])
        if fill is not None and color is None:
            plt.fill_between(range(len(list2D[thisTrial])), fill[thisTrial][0], fill[thisTrial][1], alpha = 0.2, color=color)
    
    if legend is not None:
        plt.legend(legend)
        
    plt.savefig(filename)
    if show == True:
        plt.show(block=False)
    else:
        plt.close()
    return None

def setHeaders(dataSet, prefix):
    """loop through a 1D or 2D list to generate names corresponding to a given prefix + the index of each element in the innermost list
    arguments:
    dataSet -- a 1D or 2D list to serve as an index counter. In either case, it will take the number of items within the smaller level of list
    prefix -- the base name to give to all the names that will be generated
    """
    nbOfLevels = nbOfDimensions(dataSet)
    
    headers = []
    
    if nbOfLevels == 2:
        for i in range(1,len(dataSet[0])+1):
            headers += [f'{prefix}{str(i)}']
            
    elif nbOfLevels == 1:
        for i in range(1,len(dataSet)+1):
            headers += [f'{prefix}{str(i)}']
    return headers

def msToFrames(samplingRate):
    """ give the number of frames equivalent to 1 ms"""
    ratio = 1000/samplingRate
    return ratio

def countNan(dataSet):
    """ count the number of nan in a 1D or 2D list. If a 1D list is inputed, return a single value, else, return a list"""
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 2:
        nanCounter = []
        for thisTrial in dataSet:
            thisTrialNan = 0
            for thisFrame in thisTrial:
                if str(thisFrame) == 'nan':
                    thisTrialNan += 1
            nanCounter += [thisTrialNan]
    
    elif nbOfLevels == 1:
        nanCounter = 0
        for i in dataSet:
            if str(i) == 'nan':
                nanCounter += 1
    return nanCounter

def inOut(dataSet, infLimit, supLimit):
    """ determine whether all values of a 1D or 2D list lie within the limits of given values """
    
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 1:
        for thisValue in dataSet:
            if str(thisValue) != 'nan':
                if thisValue > supLimit or thisValue < infLimit:
                    outputList = 'outOfLimits'
                    break
                else:
                    outputList = 'withinLimits'
    
    elif nbOfLevels == 2:
        outputList = []
        for thisTrial in dataSet:
            inclusionMessage = 'withinLimits'
            for thisValue in thisTrial:
                if str(thisValue) != 'nan':
                    if thisValue > supLimit or thisValue < infLimit:
                        inclusionMessage = 'outOfLimits'
                        break
            outputList += [inclusionMessage]
    
    return outputList

def pixelsToDegrees(dataSet, screenSizeCm, screenSizePx, eyeScreenDistance):
    """
    dataSet: a single value, 1D or 2D list to be converted in degrees
    screenSizeCm -- the size of the screen (in the relevant dimension) in cm
    screenSizePx -- the size of the screen (in the relevant dimension) in px
    eyeScreenDistance -- the distance separating subject's eye and the screen in cm
    """
    degPerPx = degrees(atan2(screenSizeCm/2, eyeScreenDistance)) / (screenSizePx/2)
    
    if type(dataSet) is not list:
        valueInDegrees = dataSet*degPerPx
    else:
        nbOfLevels = nbOfDimensions(dataSet)
        valueInDegrees = []
        if nbOfLevels == 1:
            for thisValue in dataSet:
                valueInDegrees += [thisValue*degPerPx]
        elif nbOfLevels == 2:
            for thisTrial in dataSet:
                thisTrialList = []
                for thisValue in thisTrial:
                    thisTrialList += [thisValue*degPerPx]
                valueInDegrees += [thisTrialList]
    return valueInDegrees

def degreesToPixels(dataSet, screenSizeCm, screenSizePx, eyeScreenDistance):
    """
    dataSet: a single value, 1D or 2D list to be converted in degrees
    screenSizeCm -- the size of the screen (in the relevant dimension) in cm
    screenSizePx -- the size of the screen (in the relevant dimension) in px
    eyeScreenDistance -- the distance separating subject's eye and the screen in cm
    """
    degPerPx = degrees(atan2(0.5*screenSizeCm, eyeScreenDistance)) / (0.5*screenSizePx)
    
    if type(dataSet) is not list:
        valueInPixels = dataSet/degPerPx
    else:
        nbOfLevels = nbOfDimensions(dataSet)
        valueInPixels = []
        if nbOfLevels == 1:
            for thisValue in dataSet:
                valueInPixels += [thisValue/degPerPx]
        elif nbOfLevels == 2:
            for thisTrial in dataSet:
                thisTrialList = []
                for thisValue in thisTrial:
                    thisTrialList += [thisValue/degPerPx]
                valueInPixels += [thisTrialList]
    return valueInPixels

def pixelsToCm(dataSet, screenSizeCm, screenSizePx):
    """
    dataSet: a single value, 1D or 2D list to be converted in degrees
    screenSizeCm -- the size of the screen (in the relevant dimension) in cm
    screenSizePx -- the size of the screen (in the relevant dimension) in px
    """
    pxPerCm = screenSizePx/screenSizeCm
    
    if type(dataSet) is not list:
        valueInCm = dataSet/pxPerCm
    
    else:
        nbOfLevels = nbOfDimensions(dataSet)
        valueInCm = []
        if nbOfLevels == 1:
            for thisValue in dataSet:
                valueInCm += [thisValue/pxPerCm]
        elif nbOfLevels == 2:
            for thisTrial in dataSet:
                thisTrialList = []
                for thisValue in thisTrial:
                    thisTrialList += [thisValue/pxPerCm]
                valueInCm += [thisTrialList]
    return valueInCm    
    
def cmToPixels(dataSet, screenSizeCm, screenSizePx):
    """
    dataSet: a single value, 1D or 2D list to be converted in degrees
    screenSizeCm -- the size of the screen (in the relevant dimension) in cm
    screenSizePx -- the size of the screen (in the relevant dimension) in px
    """
    pxPerCm = screenSizePx/screenSizeCm
    
    if type(dataSet) is not list:
        valueInPx = dataSet*pxPerCm
    
    else:
        nbOfLevels = nbOfDimensions(dataSet)
        valueInPx = []
        if nbOfLevels == 1:
            for thisValue in dataSet:
                valueInPx += [thisValue*pxPerCm]
        elif nbOfLevels == 2:
            for thisTrial in dataSet:
                thisTrialList = []
                for thisValue in thisTrial:
                    thisTrialList += [thisValue*pxPerCm]
                valueInPx += [thisTrialList]
    return valueInPx

def arbitraryToMillimeters(dataSet, scalingFactor, mode='area'):
    """ convert pupil size from a.u. to millimeters based on a scaling factor
    arguments:
    dataSet --  an int, 1d or 2d list containing the pupil size value in mm
    scaling factor -- the ratio converter from a.u. to mm, can be found using computeScalingFactor()
    mode -- the variable the arbitrary units refer to (diameter, area). default=area
    """
    if type(dataSet) is not list:
        if mode == 'diameter':
            inMillimeters = dataSet*scalingFactor
        elif mode == 'area':
            inMillimeters = sqrt(dataSet)*scalingFactor
    else:
        nbOfLevels = nbOfDimensions(dataSet)
        inMillimeters = []
        
        if nbOfLevels == 2:
            for thisTrial in dataSet:
                thisTrialConvertedValue = []
                for thisValue in thisTrial:
                    if mode == 'diameter':
                        thisTrialConvertedValue += [thisValue*scalingFactor]
                    elif mode == 'area':
                        thisTrialConvertedValue += [sqrt(thisValue)*scalingFactor]
                inMillimeters += [thisTrialConvertedValue]    
        
        elif nbOfLevels == 1:
            for thisValue in dataSet:
                if mode == 'diameter':
                    inMillimeters += [thisValue*scalingFactor]
                elif mode == 'area':
                    inMillimeters += [sqrt(thisValue)*scalingFactor]
    return inMillimeters

def computeScalingFactor(actualSize, arbitrarySize, mode='area'):
    """ computes a ratio to convert pupil size in a.u. to mm.
    argments:
    actualSize -- the actual size of the measured pupil
    arbitrarySize -- the size given for the actual size in a.u.
    mode -- the mode of pupil size recording. area, or diameter, default AREA
    """
    if mode == 'diameter':
        scalingFactor = actualSize/arbitrarySize
    elif mode == 'area':
        scalingFactor = actualSize/sqrt(arbitrarySize)
    return scalingFactor

def degreesToCm(dataSet, screenSizeCm, screenSizePx, eyeScreenDistance):
    """
    dataSet: a single value, 1D or 2D list to be converted in degrees
    screenSizeCm -- the size of the screen (in the relevant dimension) in cm
    screenSizePx -- the size of the screen (in the relevant dimension) in px
    eyeScreenDistance -- the distance separating subject's eye and the screen in cm
    """
    
    # start by converting degrees to pixels
    pixelsDataSet = degreesToPixels(dataSet, screenSizeCm, screenSizePx, eyeScreenDistance)
    
    # convert pixels in cm
    cmDataSet = pixelsToCm(pixelsDataSet, screenSizeCm, screenSizePx)
    
    return cmDataSet

def cmToDegrees(dataSet, screenSizeCm, screenSizePx, eyeScreenDistance):
    """
    dataSet: a single value, 1D or 2D list to be converted in degrees
    screenSizeCm -- the size of the screen (in the relevant dimension) in cm
    screenSizePx -- the size of the screen (in the relevant dimension) in px
    eyeScreenDistance -- the distance separating subject's eye and the screen in cm
    """
    
    # start by converting cm to pixels
    pixelsDataSet = cmToPixels(dataSet, screenSizeCm, screenSizePx)
    
    # convert pixels in degrees
    degreesDataSet = pixelsToDegrees(pixelsDataSet, screenSizeCm, screenSizePx, eyeScreenDistance)
    
    return degreesDataSet

def millimetersToArbitrary(dataSet, scalingFactor, mode='area'):
    """ convert millimeters values into arbitrary ones based on a scaling factor
    arguments:
    dataSet --  an int, 1d or 2d list containing the pupil size value in mm
    scaling factor -- the ratio converter from a.u. to mm, can be found using computeScalingFactor()
    mode -- the variable the arbitrary units refer to (diameter, area). default=area
    """
    if type(dataSet) is not list:
        if mode == 'diameter':
            inArbitrary = dataSet/scalingFactor
        elif mode == 'area':
            inArbitrary = (dataSet/scalingFactor)**2
    else:
        nbOfLevels = nbOfDimensions(dataSet)
        inArbitrary = []
        
        if nbOfLevels == 2:
            for thisTrial in dataSet:
                thisTrialConvertedValue = []
                for thisValue in thisTrial:
                    if mode == 'diameter':
                        thisTrialConvertedValue += [thisValue/scalingFactor]
                    elif mode == 'area':
                        thisTrialConvertedValue += [(thisValue/scalingFactor)**2]
                inArbitrary += [thisTrialConvertedValue]    
        
        elif nbOfLevels == 1:
            for thisValue in dataSet:
                if mode == 'diameter':
                    inArbitrary += [thisValue/scalingFactor]
                elif mode == 'area':
                    inArbitrary += [(thisValue/scalingFactor)**2]
    return inArbitrary

def transversalMean(list2D):
    """ compute the mean of sublists of floats or integers that have the same index in a 2D list
    Return a 1D list of the mean of each index"""
    allMean = []
    # loop through all columns
    for thisColumn in range(len(list2D[0])):
        thisColumnData = []
        # gather each row of the same columns
        for thisRow in range(len(list2D)):
            thisColumnData += [list2D[thisRow][thisColumn]]
        # compute the mean of each list
        allMean += [mean(thisColumnData)]
    return allMean


def distanceFromPoint(xCoord, yCoord, initialCoordinates=[0,0]):
    """ compute the absolute distance between a reference point and another using x and y coordinates
    arguments:
    xCoord -- the x coordinates of the final location (int,1d or 2d list)
    yCoord -- the y coordinates of the final location (int,1d or 2d list)
    initialCoordinates -- the reference coordinates from which a distance must be computed. default [0,0]
    """
    toSubtractX = initialCoordinates[0]
    toSubtractY = initialCoordinates[1]
    
    if type(xCoord) is not list or type(yCoord) is not list:
        xCoord -= toSubtractX
        yCoord -= toSubtractY
        hypothenuse = sqrt(xCoord**2 + yCoord**2)
    else:
        hypothenuse = []
        nbOfLevels = nbOfDimensions(xCoord)
        if nbOfLevels == 1:
            xCoord = [i-toSubtractX for i in xCoord]
            yCoord = [i-toSubtractY for i in yCoord]
            for thisValue in range(len(xCoord)):
                hypothenuse += [sqrt(xCoord[thisValue]**2 + yCoord[thisValue]**2)]
        
        elif nbOfLevels == 2:
            for thisTrial in range(len(xCoord)):
                thisTrialxCoord = []
                thisTrialyCoord = []
                thisTrialHypothenuse = []
                thisTrialxCoord += [i-toSubtractX for i in xCoord[thisTrial]]
                thisTrialyCoord += [i-toSubtractY for i in yCoord[thisTrial]]
                for thisValue in range(len(thisTrialxCoord)):
                    thisTrialHypothenuse += [sqrt(thisTrialxCoord[thisValue]**2 + thisTrialyCoord[thisValue]**2)]
                hypothenuse += [thisTrialHypothenuse]
    
    return hypothenuse

def detectSaccades(list1D, velocityThreshold=0.01, continuousFrames=1, unidirectionnal=False, saveOffsets=False, showRebounds=False):
    """ detect eye saccades based on a velocity threshold 
    arguments:
    list1D -- a 1D list of gaze coordinates
    velocityThreshold -- the velocity threshold at which a gaze displacement is considered a saccade. Default 0.01 (degrees/ms)
    continuousFrames -- the minimal number of continuous frames for which the threshold must be reached to consider that a saccade has started. 1= 2frames. Default 1
    unidirectionnal -- only detect saccades towards one direction (positive or negative values only). Default False
    saveOffsets -- output saccade offsets as a list. Default False
    showRebounds -- decide whether to keep or not saccades of opposite direction that occur just after the offset of a previous saccade (overshoot correction). Default False 
    """
    allSaccadesOnsets = []
    allSaccadesOffsets = []
    startSearch = 0
    endSearch = len(list1D) - 2
    stillSaccades = True
    
    # as long as saccades are found
    while stillSaccades == True:
        positiveOnset = None
        # look for saccades onsets and offsets
        positiveOnset = velocitySearch(list1D, range(startSearch, endSearch), velocityThreshold, continuousFrames)
        positiveOffset = velocitySearch(list1D, range(startSearch, endSearch), velocityThreshold, continuousFrames, stopAtFirst=False)
        
        if positiveOnset is None:
            stillSaccades = False
            continue
        else:
            allSaccadesOnsets += [positiveOnset]
            allSaccadesOffsets += [positiveOffset]
            startSearch = positiveOffset
    
    if unidirectionnal == False:
        stillSaccades = True
        startSearch = 0
        # as long as saccades are found
        while stillSaccades == True:
            negativeOnset = None
            # look for saccades onsets and offsets
            negativeOnset = velocitySearch(list1D, range(startSearch, endSearch), -velocityThreshold, continuousFrames)
            negativeOffset = velocitySearch(list1D, range(startSearch, endSearch), -velocityThreshold, continuousFrames, stopAtFirst=False)
            
            if negativeOnset is None:
                stillSaccades = False
                continue
            else:
                allSaccadesOnsets += [negativeOnset]
                allSaccadesOffsets += [negativeOffset]
                startSearch = negativeOffset
    
    allSaccadesOnsets.sort()
    allSaccadesOffsets.sort()
    
    if showRebounds == False:
        count=0
        while count < len(allSaccadesOnsets)-1:
            if allSaccadesOnsets[count+1] - allSaccadesOffsets[count] <= 10:
                del allSaccadesOnsets[count+1]
                del allSaccadesOffsets[count+1]
            count += 1
    
    if saveOffsets == True:
        return (allSaccadesOnsets, allSaccadesOffsets)
    else:
        return allSaccadesOnsets
    
def makeBins(dataSet, binSize):
    """ creates bins from a 1D or 2D list of int or float
    arguments:
    dataSet -- a 1D or 2D list in which bins must be made
    binSize -- the number of frames to average per bin
    """
    nbOfLevels = nbOfDimensions(dataSet)
    bins = []
    stepSize = binSize
    
    if nbOfLevels == 1:
        bins += [mean(dataSet[i:i+stepSize+1]) for i in range(0,len(dataSet), stepSize)]
        
    elif nbOfLevels == 2:
        bins += [[mean(thisList[i:i+stepSize+1]) for i in range(0,len(thisList), stepSize)] for thisList in dataSet]
    
    return bins

def interleaved(dataSet, split):
    """ divide a dataSet into 'test' and 'train' rows in an interleaved fashion. 
    Return a list of train and test sets corresponding to each interleaved possibility
    ex: data split into 3 subsets: will return 3 training sets [0,0,1],[1,0,0],[0,1,0], 1 being test and 0 training.
    arguments:
    dataSet -- a 1D or 2D list that has to be divided
    split -- the number of sets to generate
    """
    
    nbOfsplits = split
    allTrainingSets = []
    allTestSets = []
    # create the matrix of selection
    for thisCurrentTraining in range(nbOfsplits):
        selectMatrix = [0]*nbOfsplits
        selectMatrix[thisCurrentTraining] = 1
        # generate a list of the current selection list
        sequence = []
        while len(sequence) < len(dataSet):
            sequence += selectMatrix
        
        # loop into the dataset and attribute subsets
        test = []
        training = []
        for i,j in zip(sequence, dataSet):
            if i == 0:
                training += [j]
            else:
                test += [j]
        allTrainingSets += [training]
        allTestSets += [test]
    return allTrainingSets, allTestSets
