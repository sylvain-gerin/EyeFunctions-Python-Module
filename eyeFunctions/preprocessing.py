""" this module provides a series of functions used to preprocess eye tracking data recorded with Eyelink devices.It can also be used for other time series data
Dependencies: scipy.interpolate, eyeFunctions.toolBox
Author: Sylvain Gerin, sylvain.gerin@uclouvain.be 
"""

from scipy.interpolate import interp1d
from eyeFunctions.toolBox import *

def cubicSpline(list1D, onset, offset):
    """ find 4 points on a list to allow for a cubic interpolation. If one of the points is out of range, or is nan, return None
    arguments:
    list1D -- a 1D list for which points must be found to allow cubic interpolation
    onset -- the list index of the interpolation onset
    offset -- the list index of the interpolation offset
    """
    duration = offset - onset
    points = [onset-duration, onset, offset, offset+duration]
    for thisPoint in points:
        # check if those points are possible to apply
        if thisPoint < 0 or thisPoint >= len(list1D) or str(list1D[thisPoint]) == 'nan':
            points = None
    return points

def interpolate(list1D, start, stop):
    """ interpolate data of a list using a cubic interpolation if possible, and a linear one if not
    arguments:
    list1D -- the 1D list in which data must be interpolated
    start -- the list index for interpolation onset
    stop -- the list index for interpolation offset
    """
    # start by defining if a cubic interpolation is feasible
    indexPoints = cubicSpline(list1D, start, stop)
    
    # if not possible, perform a linear interpolation
    if indexPoints is None:
        # if one of the two points is nan, find the nearest outer non nan value
        if str(list1D[start]) == 'nan' and start != 0:
            for thisLowerValue in range(start,-1,-1):
                if str(list1D[thisLowerValue]) != 'nan':
                    start = thisLowerValue
                    break
        if str(list1D[stop]) == 'nan' and start != 0:
            for thisUpperValue in range(stop,len(list1D)):
                if str(list1D[thisUpperValue]) != 'nan':
                    stop = thisUpperValue
                    break
        
        indexPoints = [start, stop]
        # get the data corresponding to the points of reconstruction
        dataPoints = [list1D[i] for i in indexPoints]
        # interpolate these data
        interpolatedResult = interp1d(indexPoints, dataPoints)
    
    # if possible, perform a cubic interpolation
    else:
        # get the data corresponding to the points of reconstruction
        dataPoints = [list1D[i] for i in indexPoints]
        # interpolate these data
        interpolatedResult = interp1d(indexPoints, dataPoints, kind='cubic')
    
    # create a range of values for which interpolated data must be available
    interpolatedIndex = [i for i in range(start, stop)]
    # select the interpolated data for this range
    interpolatedData = interpolatedResult(interpolatedIndex)
    
    # put it back in the original list
    list1D[start:stop] = interpolatedData
    return None

def segmentTrials(samples, events, startMessage, endMessage, eventClockIndex=1, sampleClockIndex=0, pupilSizeIndex=3, showInfo=False, partial=False):
    """ take event and sample datasets and cluster samples based on messages in the event file
    arguments:
    samples -- the 2D list containing samples
    events -- the 2D list containing events
    startMessage -- the message from which samples are gathered
    endMessage -- the message up to which samples are gathered
    eventClockIndex -- the column index of the clock value in the event file (default=1)
    sampleClockIndex -- the column index of the clock value in the sample file (default=0)
    pupilSizeIndex -- the column index of pupil size in the sample file (default=3)
    """
    
    # retrieve all the onsets and offsets of trials in the events files
    startClockValue = []
    endClockValue = []
    # check every messages in the event file
    for messages in events:
        if partial == False:
            if startMessage in messages:
                startClockValue += [messages[eventClockIndex]] # only take the clock value of onset message
            elif endMessage in messages:
                endClockValue += [messages[eventClockIndex]] # only take the clock value of offset message
                
        elif partial == True:
            for thisMessage in messages:
                if startMessage in thisMessage:
                    startClockValue += [messages[eventClockIndex]]
                elif endMessage in thisMessage:
                    endClockValue += [messages[eventClockIndex]]
    
    # get the indexes of onset and offset messages in the sample file
    trialOnset = []
    trialOffset = []
    
    sampleClocks = getValues(samples, sampleClockIndex)
    
    for start, stop in zip(startClockValue, endClockValue):
        trialOnset += [sampleClocks.index(start)]
        trialOffset += [sampleClocks.index(stop)]
        if showInfo== True:
            print('trial', len(trialOnset))
    
    # select all samples between the limits for each trial
    allTrials = []
    for start, stop in zip(trialOnset, trialOffset):
        allTrials += [stringToFloat(samples[start:stop+1], zeroIsNan=True)]
    
    return allTrials

def downSample(dataSet, initialRate, finalRate):
    """ downsample a dataset from an initial to a final sampling rate
    arguments:
    dataSet -- the dataset to downsample, 1D or 2D list
    initialRate -- the initial sampling rate in Hz
    finalRate -- the sample rate desired in Hz for the dataset
    """
    # define the ratio at which samples will be collected
    stepSize = initialRate // finalRate
    
    nbOfLevels = nbOfDimensions(dataSet)
    downSampledData = []
    
    if nbOfLevels == 2:
        for thisTrial in dataSet:
            thisTrialdownSampledData = []
            for thisFrame in range(0,len(thisTrial),stepSize):
                thisTrialdownSampledData += [thisTrial[thisFrame]]
            downSampledData += [thisTrialdownSampledData]
    
    elif nbOfLevels == 1:
        for thisFrame in range(0,len(dataSet), stepSize):
            downSampledData += [dataSet[thisFrame]]
    
    return downSampledData

def getConditions(eventFile, independentVariable, eventMessageIndex=2):
    """ find the messages corresponding to a list of possibilities
    arguments:
    eventFile -- the dataSet in which messages must be found
    independentVariable --  a list of the messages corresponding to the levels of a given IV
    eventMessageIndex -- the column index in which messages are found in the event file
    """
    # get the informations of each trial
    conditionInfo = []
    # loop through events and take the messages of trials for each condition
    for messages in eventFile: 
        if len(messages) > eventMessageIndex: # make sure that lines without message are not taken
            if messages[eventMessageIndex] in independentVariable:
                conditionInfo += [messages[eventMessageIndex]]
    return conditionInfo

def reconstructBlinks(pupilData, xGazeData=None, yGazeData=None, velocityThreshold=3, continuousFrames=2, offsetThreshold=1, margin=10, maxDuration=500, reconstructGaze=False, countBlinks=True):
    """ interpolate pupil size data with a cubic or linear interpolation
    arguments:
    pupilData --  a 1D list of pupil size values
    xGazeData -- a 1D list of x gaze coordinates if gaze is to be reconstructed default None
    yGazeData -- a 1D list of y gaze coordinates if gaze is to be reconstructed. default None
    velocityThreshold -- the velocity threshold at which pupil size changes (both increases and decreases) indicate blinks, default 3
    continuousFrames -- the number of consecutive time points above threshold to consider that a given event occurs (blink onset, blink offset), default 2
    offsetThreshold -- the range around which velocity differences must stay to be considered as a blink offset. default 1
    margin -- the number of time points taken around the onset and offset of a blink that are considered unreliable, default 10
    maxDuration -- the maximal number of consecutive time points to interpolate, default=500
    reconstructGaze -- interpolate gaze coordinates during blinks. If not interpolated, replace with missing values. default False
    countBlinks -- decide whether to keep track of the number of detected blinks. default True
    """
    
    # this is a recursive process: it will continue until there are blinks to reconstruct
    blinkCounter =  None if countBlinks == False else 0
    
    stillBlinks = True
    while stillBlinks:
        
        # Detect blink onsets based on velocity
        thisBlinkOnset = velocitySearch(pupilData, range(len(pupilData) -2), -velocityThreshold, continuousFrames)
        # if no more blinks are found, check that no more rebounds are found
        if thisBlinkOnset is None:
            # check that no more rebounds are found, and take the first reversal
            thisReversalOnset = velocitySearch(pupilData, range(len(pupilData) -2), velocityThreshold, continuousFrames)
           
            # if no more reversals, stop the process
            if thisReversalOnset is None:
                stillBlinks = False
                continue
           
            # if a rebound is found, find an offset
            else:
                thisBlinkOffset = velocitySearch(pupilData, range(thisReversalOnset,len(pupilData) -2), offsetThreshold, continuousFrames, around=True)
                endNan = thisBlinkOffset
                
                # if no offset is found, just find the end of the reversal
                if thisBlinkOffset is None:
                    thisReversalOffset = velocitySearch(pupilData, range(thisReversalOnset,len(pupilData) -2), velocityThreshold, continuousFrames, stopAtFirst=False)
                    endNan = thisReversalOffset               
                    # if no reversal offset is found, remove all data until the end of the trial
                    if thisReversalOffset is None:
                        endNan = len(pupilData)-1
                
                # delete all data from the reversal onset to offset
                for thisValue in range(thisReversalOnset, endNan+1):
                    pupilData[thisValue] = float('nan')
                    
                    if xGazeData is not None:
                        xGazeData[thisValue] = float('nan')
                    if yGazeData is not None:
                        yGazeData[thisValue] = float('nan')
            
            stillBlinks = False
            continue
        
        # if an onset is found, look for a rebound
        if countBlinks == True:
            blinkCounter += 1
        
        # look for a rebound in the period from blink onset to the end of the recording
        startThisSearch = thisBlinkOnset
        endThisSearch = len(pupilData)-1
            
        # search rebound as the last period in which the velocity of pupil size increase exceeds the given threshold
        thisReversalOffset = velocitySearch(pupilData, range(startThisSearch, endThisSearch-1), velocityThreshold, continuousFrames, stopAtFirst=False)
        
        if thisReversalOffset == None:
            thisReversalOffset = startThisSearch
        # if no temporaly close reversal has been found, start search for a blink offset from the detected onset to avoid missing offsets
        elif(thisReversalOffset - thisBlinkOnset) > maxDuration:
            thisReversalOffset = startThisSearch
            previousReversalOffset = thisReversalOffset
            # make sure that the algorithm doesn't get stuck by trying to reconstruc always the same blink
            if thisReversalOffset == previousReversalOffset:
                thisReversalOffset = velocitySearch(pupilData, range(startThisSearch, endThisSearch-1), velocityThreshold, continuousFrames, stopAtFirst=False)
        
        # detect blink offset
        thisBlinkOffset = velocitySearch(pupilData, range(thisReversalOffset, endThisSearch-1), offsetThreshold, continuousFrames, around=True)
        # if no offset is found, assume it is because it happens at the end of recording time
        if thisBlinkOffset is None: 
            # try to check if there is a shorter continuous flat period
            thisBlinkOffset = velocitySearch(pupilData, range(thisReversalOffset, endThisSearch-1), offsetThreshold, continuousFrames//2, around=True)
            if thisBlinkOffset is None: # if it doesn't work, find the last available value or delete the end of trial
                
            # if the last value is missing, remove all values from onset to the end of recording
                if str(pupilData[len(pupilData)-1]) == 'nan':
                    for thisValue in range(thisBlinkOnset,len(pupilData)):
                        pupilData[thisValue] = float('nan')
                    # skip the interpolation
                    continue
                # if the last value is available, interpolate with that value
                else:
                    thisBlinkOffset = len(pupilData)-1
                
        # define the period to reconstruct from onset and offset
        startReconstruct = thisBlinkOnset - margin if  thisBlinkOnset - margin > 0 else 0
        endReconstruct = thisBlinkOffset + margin if thisBlinkOffset + margin < len(pupilData) else len(pupilData) -1
        
        # if the loss of signal is too long, supress the data
        if endReconstruct - startReconstruct > maxDuration:
            for thisValue in range(startReconstruct, endReconstruct+1):
                pupilData[thisValue] = float('nan')
                    
                if xGazeData is not None:
                    xGazeData[thisValue] = float('nan')
                if yGazeData is not None:
                    yGazeData[thisValue] = float('nan')
        else:
            #interpolate pupil data, and potentially gaze coordinates
            interpolate(pupilData, startReconstruct, endReconstruct)
            if reconstructGaze == True:
                if xGazeData is not None:
                    interpolate(xGazeData, startReconstruct, endReconstruct)
                if yGazeData is not None:
                    interpolate(yGazeData, startReconstruct, endReconstruct)
                
            elif reconstructGaze == False:
                for thisValue in range(startReconstruct, endReconstruct+1):
                    if xGazeData is not None:
                        xGazeData[thisValue] = float('nan')
                    if yGazeData is not None:
                        yGazeData[thisValue] = float('nan')
    
    if xGazeData is None and yGazeData is None:
        if countBlinks is False:
            return pupilData
        else:
            return pupilData, blinkCounter
    else:
        if countBlinks is False:
            return pupilData, xGazeData, yGazeData
        else:
            return pupilData, xGazeData, yGazeData, blinkCounter
    

def reconstructNan(pupilData, xGazeData=None, yGazeData=None, maxDuration=500, margin=10, reconstructGaze=False):
    """ interpolate missing data with a cubic or linear interpolation
    arguments:
    pupilData --  a 1D list of pupil size values
    xGazeData -- a 1D list of x gaze coordinates if gaze is to be reconstructed default None
    yGazeData -- a 1D list of y gaze coordinates if gaze is to be reconstructed. default None
    margin -- the number of frames taken around the onset and offset of a missing period that are considered unreliable, default 10
    maxDuration -- the maximal number of missing time points to interpolate, default 500
    reconstructGaze -- interpolate gaze coordinates. If not interpolated, keeps missing values. default False
    """
    
    stillReconstruct = True
    
    while stillReconstruct == True:
        continuousNan = 0
        nanOffset = None
        nanOnset = None
        # find continuous periods of missing values
        for thisValue in range(len(pupilData)-1):
            if str(pupilData[thisValue]) == 'nan'and thisValue != 0:
                nanOnset = thisValue - (continuousNan+1)
                continuousNan += 1
                # if it is the last missing value AND data is reconstructible, set the offset of the missing period
                if str(pupilData[thisValue+1]) != 'nan' and continuousNan < maxDuration:
                    nanOffset = thisValue+1
                    break
            else:
                continuousNan = 0
        # if no offset has been found, stop searching
        if nanOffset is None:
            stillReconstruct = False
            continue
        else:
            # set limits of reconstruction
            startReconstruct = nanOnset - margin if nanOnset - margin > 0 else 0
            endReconstruct = nanOffset + margin if nanOffset + margin < len(pupilData) else len(pupilData)-1
            
            #interpolate pupil data, and potentially gaze coordinates
            interpolate(pupilData, startReconstruct, endReconstruct)
            if reconstructGaze == True:
                if xGazeData is not None:
                    interpolate(xGazeData, startReconstruct, endReconstruct)
                if yGazeData is not None:
                    interpolate(yGazeData, startReconstruct, endReconstruct)
                    
            elif reconstructGaze == False:
                for thisValue in range(startReconstruct, endReconstruct+1):
                    if xGazeData is not None:
                        xGazeData[thisValue] = float('nan')
                    if yGazeData is not None:
                        yGazeData[thisValue] = float('nan')
    
    if xGazeData is None and yGazeData is None:
        return pupilData

    else:
        return pupilData, xGazeData, yGazeData

def backInSamples(valueDataSet, sampleDataSet, valueIndex):
    """ put newly computed values back in a larger sample dataset
    arguments:
    valueDataSet -- the newly computed values (1 or 2d list)
    sampleDataSet -- the original large dataset in which values must be replaced (2 or 3D list)
    valueIndex -- the list index of the original dataset in which new values must be put
    !!! the newly computed and original lists must have the same length !!!
    """
    nbOfLevels = nbOfDimensions(valueDataSet)
    
    newDataSet = []
    
    if nbOfLevels == 2:
        for thisTrial in range(len(valueDataSet)):
            thisTrialSample = []
            for thisFrame in range(len(valueDataSet[thisTrial])):
                # replace the original value with the new one
                sampleDataSet[thisTrial][thisFrame][valueIndex] = valueDataSet[thisTrial][thisFrame]
                thisTrialSample += [sampleDataSet[thisTrial][thisFrame]]
            newDataSet += [thisTrialSample]
    
    elif nbOfLevels == 1:
        for thisFrame in range(len(valueDataSet)):
            # replace the original value with the new one
            sampleDataSet[thisFrame][valueIndex] = valueDataSet[thisFrame]
            newDataSet += [sampleDataSet[thisFrame]]
            
    return(newDataSet)

def getBaseline(dataSet, duration, retrospective = True):
    """ compute a baseline by taking the average of a given number of time points from the beginning or the end of a list
    arguments:
    dataSet -- a  1D or 2D list for which baseline values have to be computed
    duration -- the number of time points for which the mean will be computed
    retrospective -- determine whether baseline is computed on the first frames (false) or from the last frames (true) of the list. default True
    """
    nbOfLevels = nbOfDimensions(dataSet)
    
    # convert from ms to frames
    nbOfFrames = duration
    
    if nbOfLevels == 2:
        baselineValues = []
        if retrospective == False:
            baselineValues += [mean(dataSet[thisTrial][0:nbOfFrames]) for thisTrial in range(len(dataSet))]
        elif retrospective == True:
            baselineValues += [mean(dataSet[thisTrial][len(dataSet[thisTrial])-nbOfFrames:]) for thisTrial in range(len(dataSet))]
    
    elif nbOfLevels == 1:
        baselineValues = mean(dataSet[0:nbOfFrames]) if retrospective == False else mean(dataSet[len(dataSet)-nbOfFrames:])
    
    return baselineValues

def getEventIndex(samples, events, eventMessage, eventClockIndex=1, sampleClockIndex=0):
    """ search for the clock value of a certain event, and finds the index corresponding to it in a sample file
    arguments:
    samples -- a sample file
    events -- an event file
    eventMessage -- the message to look for in the events
    eventClockIndex -- the column index for the clock value in the event file, default 1
    sampleClockIndex -- the column index for the clock value in the sample file, default 0
    """
    # retrieve clock values for each event of interest
    eventClockValue = []
    for messages in events:
        if eventMessage in messages:
            eventClockValue += [messages[eventClockIndex]]
    
    # search these clock values in the samples and get their index
    sampleClocks = getValues(samples, sampleClockIndex)
    eventIndex = [sampleClocks.index(thisClock) for thisClock in eventClockValue]
    
    return eventIndex

def smoothe(dataSet, window, nbOfParameters):
    """ smoothe data using a Savitzky-Golay filter
    arguments:
    dataSet -- the dataset to smoothe, 1D or 2D list
    window -- the number of time points in which a same filter is applied
    nbOfParameters -- the number of parameters used by the filter to set a fit. 1 linear, 2 quadratic, etc.
    """
    nbOfLevels = nbOfDimensions(dataSet)
    smoothed = []
    
    if nbOfLevels == 2:
        for thisTrial in dataSet:
            thisSmoothedTrial = []
            if all(str(values) == 'nan' for values in thisTrial):
                thisSmoothedTrial = thisTrial
            else:
                thisTrialWindow = window
                thisTrialParameters = nbOfParameters
                try:
                    thisSmoothedTrial = list(savgol_filter(thisTrial, thisTrialWindow, thisTrialParameters))
                except:
                    # remove missing data and smoothe valid data
                    validData = []
                    validFrames = []
                    validToInsert = 0
                    for thisFrame in range(len(thisTrial)):
                        if str(thisTrial[thisFrame]) != 'nan':
                            validFrames += [thisFrame]
                            validData += [thisTrial[thisFrame]]
                    # reconstruct valid data
                    try:
                        smoothedPart = list(savgol_filter(validData, thisTrialWindow, thisTrialParameters))
                    except:
                        thisTrialWindow = len(validData)
                        thisTrialParameters = thisTrialWindow - 1
                        smoothedPart = list(savgol_filter(validData, thisTrialWindow, thisTrialParameters))
                    # put the interpolated data back in the dataSet
                    for thisFrame in range(len(thisTrial)):
                        if validToInsert < len(validFrames) and thisFrame == validFrames[validToInsert]:
                            thisSmoothedTrial += [smoothedPart[validToInsert]]
                            validToInsert += 1
                        else:
                            thisSmoothedTrial += [float('nan')]
            smoothed += [thisSmoothedTrial]

    elif nbOfLevels == 1:
        if all(str(values) == 'nan' for values in dataSet):
            smoothed = dataSet
        else:
            try:
                smoothed = list(savgol_filter(dataSet, window, nbOfParameters))
            except:
                # remove missing data and smoothe valid data
                validData = []
                validFrames = []
                validToInsert = 0
                for thisFrame in range(len(dataSet)):
                    if str(dataSet[thisFrame]) != 'nan':
                        validFrames += [thisFrame]
                        validData += [dataSet[thisFrame]]
                # reconstruct valid data
                try:
                    smoothedPart = list(savgol_filter(validData, window, nbOfParameters))
                except:
                    window = len(validData)
                    nbOfParameters = window - 1
                    smoothedPart = list(savgol_filter(validData, window, nbOfParameters))
                # put the interpolated data back in the dataSet
                for thisFrame in range(len(dataSet)):
                    if validToInsert < len(validFrames) and thisFrame == validFrames[validToInsert]:
                        smoothed += [smoothedPart[validToInsert]]
                        validToInsert += 1
                    else:
                        smoothed += [float('nan')]
        
    return smoothed
