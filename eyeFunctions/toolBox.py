""" this module contains a series of functions that are useful when handling with eye-tracking data, but also when handling with lists in general.

If you encounter bugs of have any questions regarding this module, you can write an email to Sylvain Gerin: sylvain.gerin@uclouvain.be.

dependencies: matplotlib.pyplot, scipy.signal, math, csv
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

def loadFile(fileName, delimiter='\t'):
    with open(fileName, newline ='') as file:
        data = csv.reader(file, delimiter=delimiter)
        data = list(data)
    return data

def mean(list1D):
    """ Compute the mean of a list of values."""
    total = 0.0
    for i in list1D:
        if str(i) != 'nan':
            total += i
    meanValue = total/len(list1D)
    return(meanValue)

def median(list1D):
    """ compute the median of a list of values"""
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
    """ Compute the variance of a list of values."""
    meanValue = mean(list1D)
    squaredDiff = 0.0
    for i in list1D:
        squaredDiff += (i - meanValue)**2
    variance = squaredDiff / (len(list1D)-1)
    return(variance)

def getMAD(list1D, coefficient=1.4826):
    """ compute the median absolute deviation of a list of values"""
    list1D = [i for i in list1D if str(i) != 'nan']
    toSubtract = median(list1D)
    absoluteDeviations = [abs(i - toSubtract) for i in list1D]
    MAD = median(absoluteDeviations)
    MAD *=coefficient
    return MAD


def sqrt(x):
    """ Return the square root of a given number."""
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
    """ convert a 1D or 2D list of values into float. If not possible, replace the value with nan"""
    
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

def isOutlier(list1D, maxSD):
    """ get the index of trials for which a baseline pupil side exceeds a certain number of standard deviations
    arguments:
    list1D -- the 1D list of all values to consider
    maxSD -- the maximal number of standrad deviations before considering values as outliers
    """
    meanBaseline = mean(list1D)
    sdBaseline = sqrt(variance(list1D))
    upLimit = meanBaseline + maxSD*sdBaseline
    downLimit = meanBaseline - maxSD*sdBaseline
    outlierTrial = []
    
    for thisTrial in range(len(list1D)):
        if list1D[thisTrial] > upLimit or list1D[thisTrial] < downLimit:
            outlierTrial += ['outlier']
        else:
            outlierTrial += ['notOutlier']
    return(outlierTrial)

def saveList(dataSet, fileName, charSep, blockSep='\n'):
    """ Save elements of a list of 1 or 2 dimensions in a specified file.
    arguments:
    dataSet -- the v-list to save
    fileName -- the name given to the file in which the list is written
    charSep -- the separator between the list elements
    blockSep -- the separator between several blocks written on the same file (default '\n')
    """
    # Open a file in 'append' or 'write' mode
    output = open(fileName,'w')
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
                output.writelines(str(i) + charSep)
            elif i >= (len(dataSet)) - 1:
                output.writelines(str(dataSet[i]))
            output.writelines('\n')
        output.writelines(blockSep)
    else:
        print('invalid number of Dimensions of the input list')
    # Close file
    output.close()
    return None

def centerFirstFrame(dataSet):
    """subtract values to be expressed with respect to the first frame. works in 1D and 2D lists"""
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

def plotTrials(list2D, filename, legend=None, fill=None, color=None):
    """ plots a value of interest for each trial"""
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
        
    plt.savefig(filename+'.png')
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
            headers += [prefix+str(i)]
            
    elif nbOfLevels == 1:
        for i in range(1,len(dataSet)+1):
            headers += [prefix+str(i)]
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
    """ determine whether all values of a 1D or 2D list lie within the limits of given values"""
    
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 1:
        for thisValue in dataSet:
            if str(thisValue) != 'nan':
                if thisValue > supLimit or thisValue < infLimit:
                    outputList = 'outOfLimits'
                    break
                else:
                    outputList = 'in'
    
    elif nbOfLevels == 2:
        outputList = []
        for thisTrial in dataSet:
            inclusionMessage = 'in'
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
    """ compute the mean of rows in a 2D list """
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

def cluster(list2D, col, cond):
    """ Group all the rows of a 2D list for which a given column has the same value and return a list of rows (2D list).
    arguments:
    list2D -- the 2D list to search in
    col -- the column in which the condition applies
    cond -- the clustering criterion
    """
    # Create an empty list
    clusteredList = []
    # Fill it with all the rows responding to the condition
    for i in range(len(list2D)):
        if str(list2D[i][col]) == str(cond): # THE INNER LIST MUST BE STRINGS
            clusteredList += [list2D[i]]
        else:
            continue
    # Return the filled list
    return clusteredList

def distanceFromPoint(xCoord, yCoord, initialCoordinates=[0,0]):
    """ compute the absolte distance between a reference point and another using x and y coordinates
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

def detectSaccades(dataSet, velocityThreshold=0.01, continuousFrames=1, unidirectionnal=False, saveOffsets=False, showRebounds=False):
    """ detect eye saccades based on a velocity threshold 
    arguments:
    dataSet -- a 1D list of gaze coordinates
    velocityThreshold -- the velocity threshold at which a gaze displacement is considered a saccade. default 0.01 (degrees/ms)
    continuousFrames -- the minimal duration for which the threshold must be reached to consider a saccade has started
    unidirectionnal -- if you are looking only for saccades to one direction. default False
    saveOffsets -- output saccade offsets as a list
    showRebounds -- decide whether to keep or not saccades of opposite direction that occur just after the offset of a previous saccade
    """
    allSaccadesOnsets = []
    allSaccadesOffsets = []
    startSearch = 0
    endSearch = len(dataSet) - 2
    stillSaccades = True
    
    # as long as saccades are found
    while stillSaccades == True:
        positiveOnset = None
        # look for saccades onsets and offsets
        positiveOnset = velocitySearch(dataSet, range(startSearch, endSearch), velocityThreshold, continuousFrames)
        positiveOffset = velocitySearch(dataSet, range(startSearch, endSearch), velocityThreshold, continuousFrames, stopAtFirst=False)
        
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
            negativeOnset = velocitySearch(dataSet, range(startSearch, endSearch), -velocityThreshold, continuousFrames)
            negativeOffset = velocitySearch(dataSet, range(startSearch, endSearch), -velocityThreshold, continuousFrames, stopAtFirst=False)
            
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

def getValues(dataSet, variableOfInterest):
    """ gather the values of interest from the sample file (either clustered per trials or not) and return them in the same structure as the input file
    arguments:
    dataSet -- the samples in which the values are taken, can be a 1D or 2D list
    variableOfInterest -- the index of the needed value in the sample file
    """
    nbOfLevels = nbOfDimensions(dataSet)
    
    if nbOfLevels == 3:
        getSamples = []
        for thisTrial in dataSet:
            thisTrialSamples = []
            for thisFrame in thisTrial:
                thisTrialSamples += [thisFrame[variableOfInterest]]
            getSamples += [thisTrialSamples]
    
    elif nbOfLevels == 2:
        getSamples = []
        for thisFrame in dataSet:
            getSamples += [thisFrame[variableOfInterest]]
    return(getSamples)

def selectPeriod(dataSet, startPeriod, endPeriod):
    """ select only the samples from the period of interest (the recording time within a trial usually exceeds the interest period)
    arguments:
    dataSet -- the data in which trials need be shortened, a 1D or 2D list
    startPeriod -- the first frame of interest from the trial start message
    endPeriod -- the last frame of interest from the trial start message
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


def velocitySearch(list1D, rangeOfSearch, eventVelocityThreshold, eventContinuousFrames, stopAtFirst=True, around=False):
    """ iterate in a given range of data to find an event based on derivative of values, and return the index at which an event occured.
    If no event has been found, return None
    arguments:
    list1D -- a 1D list to iterate in
    rangeOfSearch -- range(onset, offset) in which an event is looked for
    eventVelocityThreshold -- the minimal difference of values between two consecutive iterations. Negative threshold indicate a decrease
    eventContinuousFrames -- the number of continuous iterations above threshold necessary to consider an event. 1=2 consecutive iterations
    stopAtFirst -- if True, stops and returns the first event found, if False, returns the last event found, default=True
    around -- if velocity must stay within a range of +- the given threshold. If True, make sure the threshold is a positive value. default=False
    """
    counter = 0
    thisEvent = None
    
    # loop through all pupil data of a given trial
    for thisFrame in rangeOfSearch:
        # check if velocity of is superior to a given threshold
        # if the event is a decrease, check for values below the threshold
        if eventVelocityThreshold < 0:
            if list1D[thisFrame+1] - list1D[thisFrame] <= eventVelocityThreshold:
                counter += 1
                # check if velocity has been above threshold for the wanted number of consecutive frames
                if counter >= eventContinuousFrames:
                    if stopAtFirst == True:
                        thisEvent = thisFrame - (counter-1)
                        break
                    else:
                        # in the case of multiple blinks, need to look for the first reversal
                        # if the event is found and will be lost at the next iteration, take this event and stop search
                        if list1D[thisFrame+2] - list1D[thisFrame+1] > eventVelocityThreshold:
                            thisEvent = thisFrame +1
                            break
                        else:
                            thisEvent = thisFrame +1
            # if the velocity is below threshold, counter is reset
            else:
                counter = 0
        elif eventVelocityThreshold > 0:
            if around == False:
                if list1D[thisFrame+1] - list1D[thisFrame] >= eventVelocityThreshold:
                    counter += 1
                    # check if velocity has been above threshold for the wanted number of consecutive frames
                    if counter >= eventContinuousFrames:
                        if stopAtFirst == True:
                            thisEvent = thisFrame - (counter-1)
                            break
                        else:
                            # in the case of multiple blinks, need to look for the first reversal
                            # if the event is found and will be lost at the next iteration, take this event and stop search
                            if list1D[thisFrame+2] - list1D[thisFrame+1] < eventVelocityThreshold:
                                thisEvent = thisFrame +1
                                break
                            else:
                                thisEvent = thisFrame +1
                # if the velocity is below threshold, counter is reset
                else:
                    counter = 0
            elif around == True:
                if abs(list1D[thisFrame+1] - list1D[thisFrame]) <= eventVelocityThreshold:
                    counter += 1
                    # check if velocity has been above threshold for the wanted number of consecutive frames
                    if counter >= eventContinuousFrames:
                        if stopAtFirst == True:
                            thisEvent = thisFrame - (counter-1)
                            break
                        else:
                            # in the case of multiple blinks, need to look for the first reversal
                            # if the event is found and will be lost at the next iteration, take this event and stop search
                            if abs(list1D[thisFrame+2] - list1D[thisFrame+1]) > eventVelocityThreshold:
                                thisEvent = thisFrame +1
                                break
                            else:
                                thisEvent = thisFrame +1
                # if the velocity is below threshold, counter is reset
                else:
                    counter = 0
    return thisEvent

def createDataBase(dataSet, nbOfTrials):
    """ combine the information of different 1D and 2D lists """
    
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
    """ compute the mean, SD and SE of a 1D or 2D list"""
    
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
    list1 -- the list containing the first operandes
    list2 -- the list containing the second operandes
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
        print('lists must be the same length')
    return finalList
