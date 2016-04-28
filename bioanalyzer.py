#!/usr/local/bin/python3
#Filename: bioanalyzer.py
#Author: Andrew Avila (Andrew.Avila@pnnl.gov)
#Description: This library provides functions to analyze the exported raw data from the Labchip GX II Bioanalyzer.

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate, stats, optimize, integrate, spatial
from functools import reduce

def GetData(rawLines, offset = 2):
    return rawLines[list(i for i in range(len(rawLines)) if rawLines[i] == "DATA")[0]+offset:]

def GetNumData(data):
    return np.array(list(float(i.split(',')[0]) for i in data)), np.array(list(float(i.split(',')[1]) for i in data))

def InterpolateData(x, y):
    return interpolate.InterpolatedUnivariateSpline(x, y)

def FindPeaks(tck, xmin, xmax, xstep = 0.01, threshold = 0):
    xnew = np.arange(xmin, xmax, xstep)
    ynew = tck(xnew)
    return np.array(list([xnew[i],ynew[i]] for i in signal.find_peaks_cwt(ynew, np.arange(xmin if xmin > 0 else 1, xmax,1)) if ynew[i] > threshold))

def GetPeaksMode(numData, modeRange = 0.20):
    return (lambda data: FindPeaks(InterpolateData(*data), int(data[0][0]), int(data[0][-1]), threshold = (lambda y: y + y*modeRange)(stats.mstats.mode(data[1], axis=None)[0][0])))(numData)

def AssignLadder(peaks, ladder = [25,50,100,150,200,300,400,500,700,850,1000]):
    return np.array(list([peaks[i][0], peaks[i][1], ladder[i]] for i in range(len(peaks)))) if len(peaks) == len(ladder) else None

def GetTopPeaks(peaks, factor = 2):
    return (lambda x: x[0:len(x)/factor])(np.array(sorted(peaks, key=lambda x: x[1], reverse = True)))

def GetSize(dataFunc, data):
    return (lambda x: np.array(list([data[i][0], data[i][1], x[i]] for i in range(len(x)))))(dataFunc(data[:,0]))

def CompareTwo(ladderLines, sampleLines):
    aLadder = AssignLadder(GetPeaksMode(GetNumData(GetData(ladderLines))))
    sampleData = GetNumData(GetData(sampleLines))
    return sampleData, GetSize(InterpolateData(aLadder[:,0], aLadder[:,2]), GetTopPeaks(GetPeaksMode(sampleData)))

def CompareRef(aLadder, sampleLines):
    return sampleData, GetSize(InterpolateData(aLadder[:,0], aLadder[:,2]), GetTopPeaks(GetPeaksMode(sampleLines)))

def AverageLadders(ladders):
    return np.array(list([np.mean(list(k[j][0] for k in ladders)), np.mean(list(k[j][1] for k in ladders)), ladders[0][j][2]] for j in range(len(ladders[0]))))

def TimeToSize(ladder, sample, values):
    return (lambda x: np.array(list([x[i], values[i]] for i in range(len(x)))))(InterpolateData(ladder[:,0], ladder[:,2])(sample))

def GetRegion(sample, lower, upper, sampley = None):
    return np.array(list(i for i in sample if i[0] > lower and i[0] < upper)) if sampley == None else np.array(list([sample[i], sampley[i]] for i in range(len(sample)) if sample[i] > lower and sample[i] < upper))

def FlatRef(sample, ref):
    return np.array(list(i-ref if i > ref else 0 for i in sample))

def SampleIntersection(a, b, c, d, start = 0, step = 0.01):
    max = a[-1] if a[-1] < c[-1] else c[-1]
    xarr = np.arange(start, max, step)
    return xarr, (lambda x, y: np.array(list(x[i] if x[i] < y[i] else y[i] for i in range(len(x)))))(InterpolateData(a, b)(xarr), InterpolateData(c, d)(xarr))

def SampleSimilarity(a, b, c, d, start = 0, end = False):
    max = (a[-1] if a[-1] < c[-1] else c[-1]) if end == False else end
    bscaled = b*np.abs(integrate.quad(InterpolateData(c, d), start, max)[0])/np.abs(integrate.quad(InterpolateData(a, b), start, max)[0])
    isect = SampleIntersection(a, bscaled, c, d)
    return np.abs(integrate.quad(InterpolateData(isect[0], isect[1]), start, max)[0])/np.abs(integrate.quad(InterpolateData(a, bscaled), start, max)[0])

def GraphSimilarity(a, b, c, d, start = 0, end = False, title = "", xaxis = "Time (s)", yaxis = "Fluorescence", sampleA = "Sample 1", sampleB = "Sample 2", ladder = []):
    max = (a[-1] if a[-1] < c[-1] else c[-1]) if end == False else end
    bscaled = b*(np.abs(integrate.quad(InterpolateData(c, d), start, max)[0])/np.abs(integrate.quad(InterpolateData(a, b), start, max)[0]))
    isect = SampleIntersection(a, bscaled, c, d)
    peaks = GetPeaksMode(isect)
    fig, ax = plt.subplots(1)
    (lambda x: ax.fill(x[:,0], x[:,1], 'r', label = sampleA, alpha = 0.3))(TimeToSize(ladder, a, bscaled)) if len(ladder) > 0 else (lambda x: ax.fill(x[:,0], x[:,1], 'r', label = sampleA, alpha = 0.3))(GetRegion(a, start, max, sampley = bscaled))
    (lambda x: ax.fill(x[:,0], x[:,1], 'b', label = sampleB, alpha = 0.3))(TimeToSize(ladder, c, d)) if len(ladder) > 0 else (lambda x: ax.fill(x[:,0], x[:,1], 'b', label = sampleB, alpha = 0.3))(GetRegion(c, start, max, sampley = d))
    (lambda x: ax.scatter(x[:,0], x[:,1], color='black', label = 'Intersection Peaks', marker = '^'))(TimeToSize(ladder, peaks[:,0], peaks[:,1])) if len(ladder) > 0 else ax.scatter(peaks[:,0], peaks[:,1], color='black', label = 'Intersection Peaks', marker = '^')
    ax.set_title(title)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.legend(loc='upper right')
    plt.show()
    return np.abs(integrate.quad(InterpolateData(isect[0], isect[1]), start, max)[0])/np.abs(integrate.quad(InterpolateData(a, bscaled), start, max)[0])

def NearestNeighbor(a, b):
    distMat = spatial.distance.cdist(a, b, 'euclidean')
    newArr = []
    for i in range(len(distMat)):
        cur = 0
        for j in range(len(distMat[i])):
            if distMat[i][j] < distMat[i][cur]:
                cur = j
        newArr += [[a[i][0], a[i][1], b[cur][0], b[cur][1], distMat[i][cur]]]
    return np.array(newArr)

def CommonNeighbors(a, b):
    return np.array(list(i for i in a if [i[2], i[3], i[0], i[1], i[4]] in b.tolist()))

def GraphNeighbors(a, b, c, d, start = 0, end = False, title = "", xaxis = "Time (s)", yaxis = "Fluorescence", sampleA = "Sample 1", sampleB = "Sample 2"):
    max = (a[-1] if a[-1] < c[-1] else c[-1]) if end == False else end
    bscaled = b*(np.abs(integrate.quad(InterpolateData(c, d), start, max)[0])/np.abs(integrate.quad(InterpolateData(a, b), start, max)[0]))
    xscaling = np.max(np.append(a, c, axis=0))
    yscaling = np.max(np.append(bscaled, d, axis=0))
    peaks_ab = GetPeaksMode([a, bscaled])
    peaks_cd = GetPeaksMode([c, d])
    neihb_ab = NearestNeighbor(np.array(list([i[0]/xscaling,i[1]/yscaling] for i in peaks_ab)), np.array(list([i[0]/xscaling,i[1]/yscaling] for i in peaks_cd)))
    neihb_cd = NearestNeighbor(np.array(list([i[0]/xscaling,i[1]/yscaling] for i in peaks_cd)), np.array(list([i[0]/xscaling,i[1]/yscaling] for i in peaks_ab)))
    common_abcd = CommonNeighbors(neihb_ab, neihb_cd)
    fig, ax = plt.subplots(1)
    (lambda x: ax.fill(x[:,0], x[:,1], 'r', label = sampleA, alpha = 0.3))(GetRegion(a, start, max, sampley = bscaled))
    (lambda x: ax.fill(x[:,0], x[:,1], 'b', label = sampleB, alpha = 0.3))(GetRegion(c, start, max, sampley = d))
    ax.scatter(peaks_ab[:,0], peaks_ab[:,1], color = 'r', marker = 'x', label = sampleA + ' Peaks')
    ax.scatter(peaks_cd[:,0], peaks_cd[:,1], color = 'b', marker = '+', label = sampleB + ' Peaks')
    list(ax.plot([i[0]*xscaling,i[2]*xscaling],[i[1]*yscaling, i[3]*yscaling], 'k') for i in common_abcd)
    ax.set_title(title)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.legend(loc='upper right')
    plt.show()
    return None

def GraphArcDis(a, b, c, d, start = 0, end = False, title = "", xaxis = "Time (s)", yaxis = "Fluorescence", sampleA = "Sample 1", sampleB = "Sample 2", ladder = []):
    max = (a[-1] if a[-1] < c[-1] else c[-1]) if end == False else end
    bscaled = b*(np.abs(integrate.quad(InterpolateData(c, d), start, max)[0])/np.abs(integrate.quad(InterpolateData(a, b), start, max)[0]))
    peaks_ab = GetPeaksMode([a, bscaled])
    peaks_cd = GetPeaksMode([c, d])
    scored = ScoreArcDisAll(np.array(list([i[0],i[1]] for i in peaks_ab)), np.array(list([i[0],i[1]] for i in peaks_cd)))
    fig, ax = plt.subplots(1)
    (lambda x: ax.fill(x[:,0], x[:,1], 'r', label = sampleA, alpha = 0.3))(TimeToSize(ladder, a, bscaled)) if len(ladder) > 0 else (lambda x: ax.fill(x[:,0], x[:,1], 'r', label = sampleA, alpha = 0.3))(GetRegion(a, start, max, sampley = bscaled))
    (lambda x: ax.fill(x[:,0], x[:,1], 'b', label = sampleB, alpha = 0.3))(TimeToSize(ladder, c, d)) if len(ladder) > 0 else (lambda x: ax.fill(x[:,0], x[:,1], 'b', label = sampleB, alpha = 0.3))(GetRegion(c, start, max, sampley = d))
    list(ax.plot([i[0],i[2]],[i[1], i[3]], color = [i[-1]**2, 1, 1 - i[-1]**2, 0.5], linewidth = (i[-1]**2)*5) for i in (lambda x, y: np.array(list([x[i][0], scored[i][1], y[i][0], scored[i][3], scored[i][4], scored[i][5]] for i in range(len(scored)))))(TimeToSize(ladder, scored[:,0], scored[:,1]), TimeToSize(ladder, scored[:,2], scored[:,3])))
    (lambda x: ax.scatter(x[:,0], x[:,1], color = 'r', marker = 'x', label = sampleA + ' Peaks'))(TimeToSize(ladder, peaks_ab[:,0], peaks_ab[:,1]))
    (lambda x: ax.scatter(x[:,0], x[:,1], color = 'b', marker = '+', label = sampleB + ' Peaks'))(TimeToSize(ladder, peaks_cd[:,0], peaks_cd[:,1]))
    ax.set_title(title)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.legend(loc='upper right')
    plt.show()
    return scored

def ScoreArcDisAll(a, b, threshold = 0.99):
    return np.array(reduce(lambda x, y: x + [y] if y not in x and np.abs(y[-1]) > threshold else x, sorted(np.append(ScoreArcDis(a,b),ScoreArcDis(b,a,reverse = True), axis=0).tolist(), key = lambda i: np.abs(i[-1])), []))
    
def ScoreArcDis(a, b, reverse = False):
    distMat = spatial.distance.cdist(a, b, 'seuclidean')
    newArr = []
    for i in range(len(distMat)):
        cur = 0
        score = 0
        cscore = 0
        for j in range(len(distMat[i])):
            cscore = (1 - np.abs(a[i][0] - b[j][0])/np.max(np.abs(b[:,0] - a[i][0]))) * ((lambda x: (2*x if x <= 0.5 else (-2*x) + 2) if x > 0 else (2*x if x >= -0.5 else (2*x) + 2))(np.arctan2(a[i][1]-b[j][1], a[i][0]-b[j][0])/np.pi))
            if np.abs(cscore) > np.abs(score):
                score = cscore
                cur = j
        newArr += [[a[i][0], a[i][1], b[cur][0], b[cur][1], distMat[i][cur], score]] if reverse == False else [[b[cur][0], b[cur][1], a[i][0], a[i][1], distMat[i][cur], score]]
    list(print(i[-1]) for i in newArr if np.abs(i[-1]) > 1)
    return np.array(newArr)

def ExpFunc(x, a, b, c):
    return a * np.exp(-b * x) + c

def LogFunc(x, a, b):
    return a + b * np.log(x)

def Poly2Func(x, a, b, c):
    return (a * x**2) + (b * x) + c

def Poly3Func(x, a, b, c, d):
    return (a * x**3) + (b * x**2) + (c * x) + d
