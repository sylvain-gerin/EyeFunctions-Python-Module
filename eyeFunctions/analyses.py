""" this module contains a series of functions that are used for time series analyses
Dependencies: csv, random, scipy.stats, numpy, eyeFunctions.toolBox
Author: Sylvain Gerin, sylvain.gerin@uclouvain.be """

import csv
import random
from scipy.stats import percentileofscore, ttest_rel, ttest_ind, ttest_1samp
import numpy as np
from eyeFunctions.toolBox import mean

def pValue(list1D, observedStat, tail='two-sided'):
    """ Compute the p-value of an observed statistic in a given dataset.
    arguments:
    list1D -- the distribution of interest, 1D list
    observedStat -- the observed score, float or integer
    tail -- the side of the distribution to look (if 'greater', return the p-value of observedStat>distribution, if 'less', the opposite, if 'two-sided', return the p-value of a bilateral test). default 'two-sided'
    """
    # get the percentile of the observed value in the distribution
    proba = percentileofscore(list1D,observedStat)
    proba = proba.item()
        
    if tail =='greater': # look to the upper end of the distribution
        thispValue = round((1-proba/100),3)
    elif tail == 'less': # look to the lower end of the distribution
        thispValue = round((proba/100),3)
    elif tail == 'two-sided': # look to both ends (assume a symetrical distribution)
        proba = percentileofscore(list1D, observedStat)
        proba = proba.item()
        if proba < 50: # if the observed value is below the distribution mean, double the probability of observed < mean
            thispValue = 2*(proba/100)
            thispValue = round(thispValue,3)
        elif proba >= 50: # if the observed score is above, double the probability of observed > mean
            thispValue = 2*(1-proba/100)
            thispValue = round(thispValue,3)
    return thispValue

def serialTtest(cond1, cond2, kind, tail='two-sided', omitNan=True):
    """ perform t-tests of 2 conditions for which a measure is taken a repeated number of times (time series)
    arguments:
    cond1 -- a 2D list for which the 1st dimension contains the values of all subjects at 1 time of observation, and each 2nd dimension refers to a different time of observation. If one-sample, the second2DList is the population mean
    cond2 -- the condition that will be compared to the 1st list, organized as the first one 
    kind -- the type of comparison that is being made, 'paired', 'independent', 'one-sample'
    tail -- the side of the distribution to look (if 'greater', return the p-value of cond1>cond2, if 'less', the opposite, if 'two-sided', return the p-value of a bilateral test). default 'two-sided'
    omitNan -- if True, nan will be omitted, else, t-tests with nan values will return to nan. default True
    """
    # create the variables of interest
    handleNan = 'omit' if omitNan == True else 'propagate'
    
    alltValues = []
    allpValues = []
    
    # perform one t-test at every time point
    for thisFrame in range(len(cond1)):
        if kind == 'paired':
            thisFrameT,thisFrameP = ttest_rel(cond1[thisFrame],cond2[thisFrame], alternative=tail, nan_policy=handleNan)
        elif kind == 'independent':
            thisFrameT, thisFrameP = ttest_ind(cond1[thisFrame],cond2[thisFrame], alternative=tail, nan_policy=handleNan)
        elif kind == 'one-sample':
            thisFrameT, thisFrameP = ttest_1samp(cond1[thisFrame],cond2, alternative=tail, nan_policy=handleNan)
        
        alltValues += [thisFrameT.item()]
        allpValues += [thisFrameP.item()]
    
    return alltValues, allpValues

def setClusters(pValues, tValues, threshold=0.05, mergeMargin=0, minClusterSize=0, descriptives=True):
    """ merge clusters of significant time points for time series
    arguments:
    pValues -- a 1D list of p-values
    tValues -- a 1D list of t-values
    threshold -- the threshold below or equal to which a p-value is considered significant
    mergeMargin -- the  maximal number of time points separating two clusters and for which clusters are to be merged. default 0
    minClusterSize -- the minimal size (in time points) for a cluster to be considered. Default 1
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
    
    # merge clusters closer than a given number of time points
    thisOnset = 0
    # check the proximity between the offset of one cluster and the onset of the next one
    while thisOnset < len(clusterOnsets) - 1:
        # if two clusters are closer than a given threshold, merge them
        if clusterOnsets[thisOnset+1] - clusterOffsets[thisOnset] <= mergeMargin:
                del clusterOffsets[thisOnset]
                del clusterOnsets[thisOnset +1]
        else:
            thisOnset += 1
    
    # remove clusters that are too short
    thisOnset = 0
    while thisOnset <= len(clusterOnsets) - 1:
        # if the duration of a cluster is less than a given threshold, delete this cluster
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
    """ randomly permute the conditions of two distributions, indepedently for each subject (made for paired permutations)
    arguments:
    cond1 -- first condition organized as a 2D list:the 1st dimension contains the values of all subjects at 1 time of observation, and each 2nd dimension refers to a different time of observation
    cond2 -- second condition organized as the first condition
    """
    # randomly decide if a subject is permuted or not => roll the dice SEPARATELY FOR EACH SUBJECT 
    selected = []
    for thisSubject in range(len(cond1[0])):
        permute = random.random()
        if permute > 0.5: # each subject has a 50% chance to be permuted
            selected += [thisSubject]
    
    # create pseudo datasets: if a subject is permuted, their values for every time point go to the other subset
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
    """ randomly permute the conditions of two distributions (made for independent permutations)
    arguments:
    cond1 -- first condition organized as a 2D list: first condition organized as a 2D list:the 1st dimension contains the values of all subjects at 1 time of observation, and each 2nd dimension refers to a different time of observation
    cond2 -- second condition organized as the first condition
    """
    
    # merge the two conditions to be permuted in a single dataset, resulting in one column per subject (regardless of condition), one row per time point
    toPermute = []
    # toPermute is the merged dataset
    for thisFrame in range(len(cond1)):
        toPermute += [cond1[thisFrame] + cond2[thisFrame]]
        
    # randomly assign subjects to either condition
    # generate n (number of subjects in cond1) random numbers between 0 and N (total number of subject in both conditions)
    selected = random.sample(range(len(cond1[0]) + len(cond2[0])),len(cond1[0]))
    # put selected subjects in the other dataset
    pseudocond1List = []
    pseudocond2List = []
    
    currentFrame = 0
    
    while currentFrame < len(toPermute):
        toFillCond1 = []
        toFillCond2 = []
        
        for subj in range(len(toPermute[currentFrame])):
            # if the subject is permuted, all their values go to the other distribution
            if subj in selected:
                toFillCond1 += [toPermute[currentFrame][subj]]
            # if the subject is not permuted, all their values stay in the same distribution
            elif subj not in selected:
                toFillCond2 += [toPermute[currentFrame][subj]]
                        
        pseudocond1List += [toFillCond1]
        pseudocond2List += [toFillCond2]
        currentFrame += 1
    
    return(pseudocond1List, pseudocond2List)
