import csv
import random
from scipy.stats import percentileofscore, ttest_rel, ttest_ind, ttest_1samp
import numpy as np
from eyeFunctions.toolBox import mean

def pValue(list1D, observedStat, tail='two-tailed'):
    """ Compute the percentile of an observed statistic in a given dataset.
    arguments:
    list1D -- the distribution of interest, 1D list
    observedStat -- the observed score, float or 1D list
    tail -- the side of the distribution to look (if 'greater', return the p-value of observedStat>distribution, if 'less', the opposite, if 'two-sided', return the p-value of a bilateral test). default 'two-sided'
    """
    proba = percentileofscore(list1D,observedStat)
    
    if tail =='greater': # look to the upper end of the distribution
        thispValue = round((1-proba/100),3)
    elif tail == 'less': # look to the lower end of the distribution
        thispValue = round((proba/100),3)
    elif tail == 'two-tailed': # look to both ends (assume a symetrical distribution)
        proba = percentileofscore(list1D, observedStat)
        if proba < 50: # if the observed value is below the distribution mean, double the probability of observed < mean
            thispValue = 2*(proba/100)
            thispValue = round(thispValue,3)
        elif proba >= 50: # if the observed score is above, double the probability of observed > mean
            thispValue = 2*(1-proba/100)
            thispValue = round(thispValue,3)               
    return thispValue

def serialTtest(first2DList, second2DList, kind, tail):
    """ perform t-tests of 2 conditions for which a measure is taken every frame
    arguments:
    first2DList -- a 2D list for which the 1st D contains the values of all subjects at 1 time of observation, and each 2nd D refers to a different time of observation. If one-sample, second2DList is the population mean
    second2DList -- the condition that will be compared to the 1st list, organized as the first one 
    kind -- the type of comparison that is being made, 'paired', 'independent', 'one-sample'
    tail -- the side of the distribution to look (if 'greater', return the p-value of observedStat>distribution, if 'less', the opposite, if 'two-sided', return the p-value of a bilateral test). default 'two-sided'
    """
    # create the variables of interest
    allpValues = []
    alltValues = []
    
    for thisFrame in range(len(first2DList)):
        if kind == 'paired':
            thisFrameT,thisFrameP = ttest_rel(first2DList[thisFrame],second2DList[thisFrame], alternative=tail)
        elif kind == 'independent':
            thisFrameT, thisFrameP = ttest_ind(first2DList[thisFrame],second2DList[thisFrame], alternative=tail)
        elif kind == 'one-sample':
            thisFrameT, thisFrameP = ttest_1samp(first2DList[thisFrame],second2DList, alternative=tail)
        
        alltValues += [thisFrameT]
        allpValues += [thisFrameP]
    
    return alltValues, allpValues

def setClusters(pValues, tValues, threshold=0.05, mergeMargin=0, minClusterSize=0, descriptives=True):
    """ set clusters of significant values for time series 
    arguments:
    pValues -- a 1D list of p-values
    tValues -- a 1D list of t-values
    threshold -- the threshold below which a p-value is considered significant
    mergeMargin -- the  maximal number of frames separating two clusters and for which clusters are to be merged. default 0
    minClusterSize -- the minimal size for a cluster to be considered. Default 1
    descriptives -- output the mean, min and max t-values of significant clusters
    """
    clusterOnsets = []
    clusterOffsets = []
    clusters = []
    sumsOfT = []
    minT = []
    maxT = []
    meanT = []
    
    for thisFrame in range(len(pValues)-1):
        # detect onsets as the first significant frame
        if thisFrame <= len(pValues)-2 and pValues[thisFrame] > threshold and pValues[thisFrame+1] <= threshold:
            clusterOnsets += [thisFrame+1]
        # detect offsets as the last significant frame
        if thisFrame <= len(pValues)-2 and pValues[thisFrame] <= threshold and pValues[thisFrame+1] > threshold:
            clusterOffsets += [thisFrame]
        # if the last frame is  part of a significant cluster, set the offset as the last frame
        if thisFrame == len(pValues)-2 and pValues[thisFrame] <= threshold and pValues[thisFrame+1] <= threshold:
            clusterOffsets += [thisFrame+1]
        # if the last frame is a single significant frame, make sure that it is encoded as an offset
        if thisFrame == len(pValues)-2 and pValues[thisFrame] > threshold and pValues[thisFrame+1] <= threshold:
            clusterOffsets += [thisFrame+1]
    
    # merge clusters closer than a given number of frames
    thisOnset = 0
    # check the proximity between the offset of one cluster and the onset of the next one
    while thisOnset < len(clusterOnsets) - 1:
        if clusterOnsets[thisOnset+1] - clusterOffsets[thisOnset] <= mergeMargin:
                del clusterOffsets[thisOnset]
                del clusterOnsets[thisOnset +1]
        else:
            thisOnset += 1
    
    # remove clusters that are too short
    thisOnset = 0
    
    while thisOnset <= len(clusterOnsets) - 1:
        if clusterOffsets[thisOnset] - clusterOnsets[thisOnset] < minClusterSize:
            del clusterOnsets[thisOnset]
            del clusterOffsets[thisOnset]
        else:
            thisOnset += 1
    
    # gather onset and offset in a single list
    for thisCluster in range(len(clusterOnsets)):
        clusters += [[clusterOnsets[thisCluster], clusterOffsets[thisCluster]]]
    
    # if no cluster is found, return nothing
    if len(clusters) != 0:
        for thisCluster in clusters:
            sumsOfT += [sum(tValues[thisCluster[0]:thisCluster[1]+1])]
            meanT += [mean(tValues[thisCluster[0]:thisCluster[1]+1])]
            minT += [min(tValues[thisCluster[0]:thisCluster[1]+1])]
            maxT += [max(tValues[thisCluster[0]:thisCluster[1]+1])]
    else:
        clusters = sumsOfT = meanT = minT = maxT = None
    
    if descriptives == True:
        return clusters, sumsOfT, meanT, minT, maxT
    else:
        return clusters, sumsOfT

def pairedRandomAssign(cond1, cond2):
    """ randomly permute the conditions of two distributions, indepedently for each participant.
    arguments:
    cond1 -- first condition organized as a 2D list (1row/frame, 1column/subject)
    cond2 -- second condition organized as a 2D list (1row/frame, 1column/subject)
    """
    # randomly decide if a participant is permuted or not => roll the dice SEPARATELY FOR EACH PARTICIPANT 
    selected = []
    for i in range(len(cond1[0])):
        permute = random.random()
        if permute > 0.5: # each participant has a 50% chance to be permuted
            selected += [i]
    
    # create pseudo datasets: if participant is permuted, their values go to the other subset
    pseudoCond1 = []
    pseudoCond2 = []
    
    currentFrame = 0
    
    while currentFrame < len(cond1):
        toFillCond1 = []
        toFillCond2 = []
        
        for subj in range(len(cond1[currentFrame])):
            # if the subject is permuted, all their frames go to the other distribution
            if subj in selected:
                toFillCond1 += [cond2[currentFrame][subj]]
                toFillCond2 += [cond1[currentFrame][subj]]
            # if the subject is not permuted, all their frames stay in the same distribution
            else:
                toFillCond1 += [cond1[currentFrame][subj]]
                toFillCond2 += [cond2[currentFrame][subj]]
                
        pseudoCond1 += [toFillCond1]
        pseudoCond2 += [toFillCond2]
        currentFrame += 1
        
    return(pseudoCond1, pseudoCond2)

def randomAssign(cond1, cond2):
    """ randomly permute the conditions of two distributions.
    arguments:
    cond1 -- first condition organized as a 2D list (1row/frame, 1column/subject)
    cond2 -- second condition organized as a 2D list (1row/frame, 1column/subject)
    """
    
    # merge the two conditions of a same finger in a single dataset, with 1 line per subject
    toPermute = []
    # toPermute is the merged dataset
    for i in range(len(cond1)):
        toPermute += [cond1[i] + cond2[i]]
        
    # randomly assign participants to either condition
    # generate n (n for cond1) random numbers between 0 and N (full sample)
    selected = random.sample(range(len(cond1[0]) + len(cond2[0])),len(cond1[0]))
    # generate pseudo datasets: permuted data go to the 'wrong' dataset
    pseudocond1List = []
    pseudocond2List = []
    
    currentFrame = 0
    
    while currentFrame < len(toPermute):
        toFillDown = []
        toFillUp = []
        
        for subj in range(len(toPermute[currentFrame])):
            # if the subject is permuted, all their frames go to the other distribution
            if subj in selected:
                toFillDown += [toPermute[currentFrame][subj]]
            # if the subject is not permuted, all their frames stay in the same distribution
            elif subj not in selected:
                toFillUp += [toPermute[currentFrame][subj]]
                        
        pseudocond1List += [toFillDown]
        pseudocond2List += [toFillUp]
        currentFrame += 1
    
    return(pseudocond1List, pseudocond2List)
