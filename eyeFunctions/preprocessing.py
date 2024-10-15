""" this module provides a series of custom-made functions that are useful when preprocessing eye tracking data recorded with Eyelink devices.

If you encounter bugs of have any questions regarding this module, you can write an email to Sylvain Gerin: sylvain.gerin@uclouvain.be.

dependencies: scipy.interpolate, eyeFunctions.toolBox
"""

from scipy.interpolate import interp1d
from eyeFunctions.toolBox import *

def cubicSpline(list1D, onset, offset):
    """ finds 4 points on a list to allow for a cubic interpolation. If one of the points is out of range, or is nan, return None"""
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
    list1D -- the 1D list in which data will be interpolated
    start -- the onset of the range for which data will be interpolated
    stop -- the offset of the range for which data will be interpolated
    """
    # start by defining if a cubic interpolation is feasible
    indexPoints = cubicSpline(list1D, start, stop)
    
    # if not possible, linear interpolation
    if indexPoints is None:
        indexPoints = [start, stop]
        # get the data corresponding to the points of reconstruction
        dataPoints = [list1D[i] for i in indexPoints]
        # interpolate these data
        interpolatedResult = interp1d(indexPoints, dataPoints)
    
    # if possible, cubic interpolation
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

def segmentTrials(samples, events, startMessage, endMessage, eventClockIndex=1, sampleClockIndex=0, pupilSizeIndex=3, showInfo=False):
    """ take an event and a sample datasets and cluster samples per trial based on messages in the event file
    arguments:
    samples -- the list containing samples
    events -- the list containing events
    startMessage -- the message indicating trial onset
    endMessage -- the message indicating trial offset
    eventClockIndex -- the column index of the clock value in the event file (default=1)
    sampleClockIndex -- the column index of the clock value in the sample file (default=0)
    pupilSizeIndex -- the column index of pupil size in the sample file (default=3)
    """
    # retrieve all the onsets and offsets of trials in the events files
    startClockValue = []
    endClockValue = []
    for messages in events:
        if startMessage in messages:
            startClockValue += [messages[eventClockIndex]] # only take the clock value of trial onset
        elif endMessage in messages:
            endClockValue += [messages[eventClockIndex]] # only take the clock value of trial offset
    
    # loop through samples and segment each trial
    trialOnset = []
    trialOffset = []
    startNextSearch = 0
    
    for currentIdx in range(len(startClockValue)):
        thisTrialOnset = []
        thisTrialOffset = []
        for currentSample in samples[startNextSearch:]:
            if float(startClockValue[currentIdx]) == float(currentSample[sampleClockIndex]):
                thisTrialOnset += [samples.index(currentSample)]
            if float(endClockValue[currentIdx]) == float(currentSample[sampleClockIndex]):
                thisTrialOffset += [samples.index(currentSample)]
                startNextSearch = samples.index(currentSample) # made to fasten the search
                break
        trialOnset += thisTrialOnset
        trialOffset += thisTrialOffset
        if showInfo== True:
            print('trial', currentIdx)
    
    # select all samples between the limits for each trial
    print('loading samples')
    allTrials = []
    for thisTrial in range(len(trialOnset)):
        thisTrialSamples = [columns[:pupilSizeIndex+1] for columns in samples[trialOnset[thisTrial]:trialOffset[thisTrial]+1]]
        thisTrialAllSamples = []
        for thisFrame in thisTrialSamples:
            try:
                thisFrame = [float(thisValue) for thisValue in thisFrame]
            except: # if no value, is '.' or '0' in the original file
                for thisValue in range(sampleClockIndex+1,len(thisFrame)):
                    thisFrame[thisValue] = float('nan')
            thisTrialAllSamples += [thisFrame] # is the same as before, but with floats and nan instead of strings
        allTrials += [thisTrialAllSamples]
    
    return allTrials

def downSample(dataSet, initialRate, finalRate):
    """ downsample a sample dataset from an initial to a final sample rate
    arguments:
    dataSet -- the dataset to downsample
    initialRate -- the initial sampling rate in Hz (generally the one of the eye-tracker)
    finalRate -- the sample rate desired in Hz for the dataset
    """
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
    """ get the conditions of a given independent variable, for each trial
    arguments:
    eventFile -- the dataSet containing all event messages
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

def reconstructBlinks(pupilData, xGazeData, yGazeData, velocityThreshold=3, continuousFrames=3, offsetThreshold=1, margin=10, maxDuration=500, reconstructGaze=False, countBlinks=True):
    """ interpolate pupil size data with a cubic-spline or linear interpolation (advised on smoothed data)
    arguments:
    pupilData --  a 1D list of pupil size
    xGazeData -- a 1D list of x gaze coordinates
    yGazeData -- a 1D list of y gaze coordinates
    velocityThreshold -- the velocity threshold at which pupil size changes (both increases and decreases) reflect blinks, default=3
    continuousFrames -- the number of consecutive frames above threshold to consider that a given event occurs (blink onset, blink offset), default=3
    offsetThreshold -- the range around which velocity differences must stay to be considered as a blink offset. default=1
    margin -- the number of frames taken around the onset and offset of a blink that are considered unreliable, default=10
    maxDuration -- the maximal duration to interpolate in frames, default=500
    reconstructGaze -- interpolate gaze coordinates during blinks. If not interpolated, replace with missing values. default=False
    countBlinks -- decide whether to keep track of the number of blinks per trial or not. default=True
    """
    workingPupil = pupilData.copy()
    workingGazeX = xGazeData.copy()
    workingGazeY = yGazeData.copy()
    
    # this is a recursive process: it will continue until there are blinks to reconstruct
    blinkCounter=None
    if countBlinks == True:
        blinkCounter = 0
    stillBlinks = True
    while stillBlinks == True:
        
        # Detect blink onsets based on velocity
        
        thisBlinkOnset = velocitySearch(workingPupil, range(len(workingPupil) -2), -velocityThreshold, continuousFrames)
        # if no more blinks are found:
        # check that no more rebounds are found
        if thisBlinkOnset is None:
            # check that no more rebounds are found, with a more laxist continuousFrames criterion and take the first found reversal
            thisReversalOnset = velocitySearch(workingPupil, range(len(workingPupil) -2), velocityThreshold, continuousFrames)
            # if no more reversals, stop the process
            if thisReversalOnset is None:
                stillBlinks = False
                continue
            # if a rebound is found, find an offset
            else:
                thisBlinkOffset = velocitySearch(workingPupil, range(thisReversalOnset,len(workingPupil) -2), offsetThreshold, continuousFrames, around=True)
                endNan = thisBlinkOffset
                
                # if no offset is found, just find the end of the reversal
                if thisBlinkOffset is None:
                    thisReversalOffset = velocitySearch(workingPupil, range(thisReversalOnset,len(workingPupil) -2), velocityThreshold, continuousFrames, stopAtFirst=False)
                    endNan = thisReversalOffset
                    # if no reversal offset is found, remove all data until the end of the trial
                    if thisReversalOffset is None:
                        endNan = len(workingPupil)-1
                # delete all data from the reversal onset to offset
                for thisValue in range(thisReversalOnset, endNan+1):
                    workingPupil[thisValue] = float('nan')
                    workingGazeX[thisValue] = float('nan')
                    workingGazeY[thisValue] = float('nan')
            stillBlinks = False
            continue
        
        # detect rebound
        if countBlinks == True:
            blinkCounter += 1
        # set the search period
        startThisSearch = thisBlinkOnset
        endThisSearch = len(workingPupil)-1
        
        # if recording stops before the end of the search period, just look as far as data is available
        if endThisSearch > len(workingPupil) + 1:
            endThisSearch = len(workingPupil) - 2
        
        # search rebound as the last period in which the velocity of increase of pupil size exceeds the given threshold
        thisReversalOffset = velocitySearch(workingPupil, range(startThisSearch, endThisSearch-1), velocityThreshold, continuousFrames, stopAtFirst=False)
        # if no temporaly close reversal has been found, start search for a blink offset from the detected onset
        if thisReversalOffset == None or (thisReversalOffset - thisBlinkOnset) > maxDuration:
            thisReversalOffset = startThisSearch
        
        # detect blink offset
        thisBlinkOffset = velocitySearch(workingPupil, range(thisReversalOffset, endThisSearch-1), offsetThreshold, continuousFrames, around=True)
        # if no offset is found, assume it is because it happens at the end of recording time
        if thisBlinkOffset is None: 
            # try to check if there is a shorter continuous flat period
            thisBlinkOffset = velocitySearch(workingPupil, range(thisReversalOffset, endThisSearch-1), offsetThreshold, continuousFrames//2, around=True)
            if thisBlinkOffset is None: # if it doesn't work, find the last available value or delete the end of trial
                # if the last value is missing, remove all values from onset to the end of recording
                if str(workingPupil[len(workingPupil)-1]) == 'nan':
                    for thisValue in range(thisBlinkOnset,len(workingPupil)):
                        workingPupil[thisValue] = float('nan')
                    # skip the interpolation
                    continue
                # if the last value is available, interpolate from that value
                else:
                    thisBlinkOffset = len(workingPupil)-1
                
        # define the period to reconstruct from onset and offset
        
        if thisBlinkOnset - margin > 0:
            startReconstruct = thisBlinkOnset - margin
        else:
            startReconstruct = 0
        
        if thisBlinkOffset + margin < len(workingPupil):
            endReconstruct = thisBlinkOffset + margin
        else:
            endReconstruct = len(workingPupil) -1
        
        # if the loss of signal is too long, supress the data
        if endReconstruct - startReconstruct > maxDuration:
            for thisValue in range(startReconstruct, endReconstruct+1):
                workingPupil[thisValue] = float('nan')
                workingGazeX[thisValue] = float('nan')
                workingGazeY[thisValue] = float('nan')
        else:
            #interpolate pupil data, and potentially gaze coordinates
            interpolate(workingPupil, startReconstruct, endReconstruct)
            if reconstructGaze == True:
                interpolate(workingGazeX, startReconstruct, endReconstruct)
                interpolate(workingGazeY, startReconstruct, endReconstruct)
            elif reconstructGaze == False:
                for thisValue in range(startReconstruct, endReconstruct+1):
                    workingGazeX[thisValue] = float('nan')
                    workingGazeY[thisValue] = float('nan')
    
    return (workingPupil, workingGazeX, workingGazeY, blinkCounter)

def reconstructNan(pupilData, xGazeData, yGazeData, maxDuration=500, margin=10, reconstructGaze=False):
    """ interpolate missing data with a cubic-spline or linear interpolation
    arguments:
    pupilData --  a 1D list of pupil size
    xGazeData -- a 1D list of gaze x coordinates
    yGazeData -- a 1D list of gaze y coordinates
    margin -- the number of frames taken around the onset and offset of a missing period that are considered unreliable, default=10
    maxDuration -- the maximal number of frames to interpolate, default=500
    reconstructGaze -- interpolate gaze coordinates during blinks. If not interpolated, replace with missing values. default=False
    """
    workingPupil = pupilData.copy()
    workingGazeX = xGazeData.copy()
    workingGazeY = yGazeData.copy()
    
    stillReconstruct = True
    
    while stillReconstruct == True:
        continuousNan = 0
        nanOffset = None
        nanOnset = None
        # find continuous missing values
        for thisValue in range(len(workingPupil)-1):
            if str(workingPupil[thisValue]) == 'nan':
                nanOnset = thisValue - (continuousNan+1)
                continuousNan += 1
                # if it is the last missing value AND data is reconstructible, set the onset of the missing period
                if str(workingPupil[thisValue+1]) != 'nan' and continuousNan < maxDuration:
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
            if nanOnset - margin > 0:
                startReconstruct = nanOnset - margin
            else:
                startReconstruct = 0
            
            if nanOffset + margin < len(workingPupil):
                endReconstruct = nanOffset + margin
            else:
                endReconstruct = len(workingPupil)-1
                
            #interpolate pupil data, and potentially gaze coordinates
            interpolate(workingPupil, startReconstruct, endReconstruct)
            if reconstructGaze == True:
                interpolate(workingGazeX, startReconstruct, endReconstruct)
                interpolate(workingGazeY, startReconstruct, endReconstruct)
            elif reconstructGaze == False:
                for thisValue in range(startReconstruct, endReconstruct+1):
                    workingGazeX[thisValue] = float('nan')
                    workingGazeY[thisValue] = float('nan')
    
    return (workingPupil, workingGazeX, workingGazeY)

def backInSamples(valueDataSet, sampleDataSet, valueIndex):
    """ put newly computed values back in a larger sample dataset
    arguments:
    valueDataSet -- the newly computed values, separated by trials (2d list)
    sampleDataSet -- the original large dataset in which values must be replaced (3D list)
    valueIndex -- the column index of the original dataset in which new values must be put
    !!! the newly computed and original lists must have the same length !!!
    """
    nbOfLevels = nbOfDimensions(valueDataSet)
    
    thisSessionSample = []
    
    if nbOfLevels == 2:
        for thisTrial in range(len(valueDataSet)):
            thisTrialSample = []
            for thisFrame in range(len(valueDataSet[thisTrial])):
                # replace the original value with the interpolated one
                sampleDataSet[thisTrial][thisFrame][valueIndex] = valueDataSet[thisTrial][thisFrame]
                thisTrialSample += [sampleDataSet[thisTrial][thisFrame]]
            thisSessionSample += [thisTrialSample]
    elif nbOfLevels == 1:
        for thisFrame in range(len(valueDataSet)):
            # replace the original value with the interpolated one
            sampleDataSet[thisFrame][valueIndex] = valueDataSet[thisFrame]
            thisSessionSample += [sampleDataSet[thisFrame]]
            
    return(thisSessionSample)

def getBaseline(dataSet, duration, retrospective = True):
    """ get the baseline pupil size by taking the average of a given number of frames from the beginning of the trial/the end of the previous
    arguments:
    dataSet -- a  1D or 2D list for which baseline values have to be computed
    duration -- the nuber of frames for which the mean will be computed
    retrospective -- determine whether baseline is computed on the first frames (false) of each trial, or from the last frame of this trial (true). default = (True)
    """
    nbOfLevels = nbOfDimensions(dataSet)
    
    thisSessionBaseline = []
    
    # convert from ms to frames
    nbOfFrames = duration
    
    if nbOfLevels == 2:
        if retrospective == False:
            for thisTrial in range(len(dataSet)):
                thisTrialBaseline = dataSet[thisTrial][0:nbOfFrames]
                thisSessionBaseline += [mean(thisTrialBaseline)]
        elif retrospective == True:
            for thisTrial in range(len(dataSet)):
                thisTrialBaseline = dataSet[thisTrial][len(dataSet[thisTrial])-nbOfFrames:]
                thisSessionBaseline += [mean(thisTrialBaseline)]
    elif nbOfLevels == 1:
        if retrospective == False:
            thisTrialBaseline = dataSet[0:nbOfFrames]
            thisSessionBaseline += [mean(thisTrialBaseline)]
        elif retrospective == True:
            thisTrialBaseline = dataSet[len(dataSet)-nbOfFrames:]
            thisSessionBaseline += [mean(thisTrialBaseline)]
    
    return thisSessionBaseline

def getEventIndex(samples, events, eventMessage, eventClockIndex=1, sampleClockIndex=0):
    """ search for the clock value of a certain event, and finds the index corresponding to it in a sample file
    arguments:
    samples -- a sample file
    events -- an event file
    eventMessage -- the message you are looking for in the events
    eventClockIndex -- the column index for the clock value in the event file, default=1
    sampleClockIndex -- the column index for the clock value in the sample file, default=0
    """
    # retrieve clock values for each event of interest 
    eventClockValue = []
    for messages in events:
        if eventMessage in messages:
            eventClockValue += [messages[eventClockIndex]]
    
    # search these clock values in the samples and get their index
    startNextSearch = 0
    eventIndex = []
    
    for currentIdx in range(len(eventClockValue)):
        for currentSample in samples[startNextSearch:]:
            if float(eventClockValue[currentIdx] )== float(currentSample[sampleClockIndex]):
                thisEventIndex = samples.index(currentSample)
                startNextSearch = samples.index(currentSample)
                break
        eventIndex += [thisEventIndex]
    return eventIndex

def smoothe(dataSet, window, nbOfParameters):
    """ smoothe data using a Savitzky-Golay filter
    arguments:
    dataSet -- the dataset to smoothe, can be a 1D or 2D list
    window -- the number of frames in which a same filter is applied
    nbOfParameters -- the number of parameters used by the filter to set a fit. 1=linear, 2=quadratic, etc.
    """
    nbOfLevels = nbOfDimensions(dataSet)
    smoothed = []
    
    if nbOfLevels == 2:
        for thisTrial in dataSet:
            try:
                thisSmoothedTrial = list(savgol_filter(thisTrial, window, nbOfParameters))
            except:
                   # remove missing data and smoothe valid data
                   thisSmoothedTrial = []
                   validData = []
                   validFrames = []
                   validToInsert = 0
                   for thisFrame in range(len(thisTrial)):
                       if str(thisTrial[thisFrame]) != 'nan':
                           validFrames += [thisFrame]
                           validData += [thisTrial[thisFrame]]
                   # reconstruct valid data
                   smoothedPart = list(savgol_filter(validData, window, nbOfParameters))
                   # put the interpolated data back in the dataSet
                   for thisFrame in range(len(thisTrial)):
                       if validToInsert < len(validFrames) and thisFrame == validFrames[validToInsert]:
                           thisSmoothedTrial += [smoothedPart[validToInsert]]
                           validToInsert += 1
                       else:
                           thisSmoothedTrial += [float('nan')]
            smoothed += [thisSmoothedTrial]              
    
    
    elif nbOfLevels == 1:
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
            smoothedPart = list(savgol_filter(validData, window, nbOfParameters))
            # put the interpolated data back in the dataSet
            for thisFrame in range(len(dataSet)):
                if validToInsert < len(validFrames) and thisFrame == validFrames[validToInsert]:
                    smoothed += [smoothedPart[validToInsert]]
                    validToInsert += 1
                else:
                    smoothed += [float('nan')]
    
    return smoothed
