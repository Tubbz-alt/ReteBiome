#!/usr/local/bin/python3
#Filename: analysis_service.py
#Author: Andrew Avila (Andrew.Avila@pnnl.gov)
#Description: This script performs analysis procedures in a parallel manner using a database queue.

from IPython.parallel import Client, require
from time import sleep
import mb_database as mb
import sm_database as sm


def GetData(rawLines, offset = 2):
    return rawLines[list(i for i in range(len(rawLines)) if rawLines[i] == "DATA")[0]+offset:]

def GetNumData(data):
    return numpy.array(list(float(i.split(',')[0]) for i in data)), numpy.array(list(float(i.split(',')[1]) for i in data))

def InterpolateData(x, y):
    return interpolate.InterpolatedUnivariateSpline(x, y)

def SampleIntersection(a, b, c, d, start = 0, step = 0.01):
    mav = a[-1] if a[-1] < c[-1] else c[-1]
    xarr = numpy.arange(start, mav, step)
    return xarr, (lambda x, y: numpy.array(list(x[i] if x[i] < y[i] else y[i] for i in range(len(x)))))(InterpolateData(a, b)(xarr), InterpolateData(c, d)(xarr))

def GetPeaksMode(numData, modeRange = 0.20):
    return (lambda data: FindPeaks(InterpolateData(*data), int(data[0][0]), int(data[0][-1]), threshold = (lambda y: y + y*modeRange if y > 20 else 20)(stats.mstats.mode(data[1], axis=None)[0][0])))(numData)

def FindPeaks(tck, xmin, xmax, threshold = 0):
    xnew = numpy.arange(xmin, xmax, numpy.around((xmax - xmin) / 7000, 2))
    ynew = tck(xnew)
    return numpy.array(list([xnew[i],ynew[i]] for i in signal.find_peaks_cwt(ynew, numpy.arange(xmin if xmin > 0 else 1, xmax, 1)) if ynew[i] > threshold))

def AssignLadder(peaks, ladder = [25,50,100,150,200,300,400,500,700,850,1000]):
    return numpy.array(list([peaks[i][0], peaks[i][1], ladder[i]] for i in range(len(peaks[1:-1])))) if len(peaks[1:-1]) == len(ladder) else None

def AverageLadders(ladders):
    return numpy.array(list([numpy.mean(list(k[j][0] for k in ladders)), numpy.mean(list(k[j][1] for k in ladders)), ladders[0][j][2]] for j in range(len(ladders[0]))))

def TimeToSize(ladder, sample, values):
    return (lambda x: numpy.array(sorted(list([x[i], values[i]] for i in range(len(x))), key = lambda x : x[0])))(InterpolateData(ladder[:,0], ladder[:,2])(sample))

def GetRegion(sample, lower, upper, sampley = None, zero = False):
    return numpy.array(sorted(list([i[0], i[1] if i[1] > 0 else (i[1] if zero == False else 0)] for i in sample if i[0] > lower and i[0] < upper), key = lambda x: x[0])) if sampley == None else numpy.array(sorted(list([sample[i], sampley[i] if sampley[i] > 0 else (sampley[i] if zero == False else 0)] for i in range(len(sample)) if sample[i] > lower and sample[i] < upper), key = lambda x: x[0]))

def ScoreArcDisAll(a, b, threshold = 0):
    return numpy.array(reduce(lambda x, y: x + [y] if y not in x and numpy.abs(y[-1]) > threshold else x, sorted(numpy.append(ScoreArcDis(a,b),ScoreArcDis(b,a,reverse = True), axis=0).tolist(), key = lambda i: numpy.abs(i[-1])), []))
    
def ScoreArcDis(a, b, reverse = False):
    distMat = spatial.distance.cdist(a, b, 'seuclidean')
    newArr = []
    for i in range(len(distMat)):
        cur = 0
        score = 0
        cscore = 0
        for j in range(len(distMat[i])):
            cscore = (1 - numpy.abs(a[i][0] - b[j][0])/numpy.max(numpy.abs(b[:,0] - a[i][0]))) * ((lambda x: (2*x if x <= 0.5 else (-2*x) + 2) if x > 0 else (2*x if x >= -0.5 else (2*x) + 2))(numpy.arctan2(a[i][1]-b[j][1], a[i][0]-b[j][0])/numpy.pi))
            if numpy.abs(cscore) > numpy.abs(score):
                score = cscore
                cur = j
        newArr += [[a[i][0], a[i][1], b[cur][0], b[cur][1], distMat[i][cur], score]] if reverse == False else [[b[cur][0], b[cur][1], a[i][0], a[i][1], distMat[i][cur], score]]
    return numpy.array(newArr)

def DefineRegions(samp_a, samp_b, ab_pdif, nthreshold = 5, threshold = 0.5):
    stack = []
    orient = False
    final = []
    orientF = lambda x: 0 if samp_a[1][x] < nthreshold and samp_b[1][x] < nthreshold else (1 if ab_pdif[1][x] > threshold else -1)
    for i in range(len(ab_pdif[0])):
        if orient == False and len(stack) == 0:
            stack += [ab_pdif[0][i]]
            orient = orientF(i)
        else:
            if orient == orientF(i):
                stack += [ab_pdif[0][i]]
            else:
                final += [[stack[0], stack[-1], orient]]
                stack = [ab_pdif[0][i]]
                orient = orientF(i)
    final += [[stack[0], stack[-1], orient]]
    return final

def ConnectParallel(prf = 'default'):
    return Client(profile = prf)

def wrapprint(i):
    global wrapvar
    wrapvar = i
    print(i)
    return i

@require(GetData, GetNumData, InterpolateData, SampleIntersection, GetPeaksMode, ScoreArcDisAll, ScoreArcDis, FindPeaks, TimeToSize, AverageLadders, AssignLadder)
def GraphArcDis(Samples, Primers_ID_F, Primers_ID_R):
    svars = list(list(json.loads(j[3]) for j in sm.SampleVariables(dbs, sampleID = i) if (lambda x: True if x["Primers_ID_F"] == Primers_ID_F and x["Primers_ID_R"] == Primers_ID_R else False)(json.loads(j[3])))[0] for i in Samples)
    odata = list(GetNumData(GetData(open(svars[i]['Path'] + svars[i]['File_Prefix'] + svars[i]['Data_RC'] + svars[i]['File_Suffix']).read().splitlines())) for i in range(len(Samples)))
    mval = max(list(max(i[0]) for i in odata))
    sdata = [odata[0]] + list([i[0], i[1]*(numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), 0, mval)[0])/numpy.abs(integrate.quad(InterpolateData(i[0], i[1]), 0, mval)[0]))] for i in odata[1:])
    peaks = list(GetPeaksMode([sdata[i][0], sdata[i][1]]) for i in range(len(sdata)))
    print(peaks)
    scored = (lambda x, y: {json.dumps(y[i]) : ScoreArcDisAll(x[i][0], x[i][1]).tolist() for i in range(len(x))})(list(combinations(peaks, 2)), list(combinations(Samples, 2)))
    return json.dumps(scored)
 
@require(GetData, GetNumData, InterpolateData, SampleIntersection)
def GraphSimilarity(Samples, Primers_ID_F, Primers_ID_R):
    svars = list(list(json.loads(j[3]) for j in sm.SampleVariables(dbs, sampleID = i, ) if (lambda x: True if "Primers_ID_F" in x.keys() and x["Primers_ID_F"] == Primers_ID_F and x["Primers_ID_R"] == Primers_ID_R else False)(json.loads(j[3])))[0] for i in Samples)
    odata = list(GetNumData(GetData(open(svars[i]['Path'] + svars[i]['File_Prefix'] + svars[i]['Data_RC'] + svars[i]['File_Suffix']).read().splitlines())) for i in range(len(Samples)))
    mval = max(list(max(i[0]) for i in odata))
    sdata = list([i[0], i[1]*(numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), 0, mval)[0])/numpy.abs(integrate.quad(InterpolateData(i[0], i[1]), 0, mval)[0]))] for i in odata[1:])
    isect = reduce(lambda x, y: SampleIntersection(x[0], x[1], y[0], y[1]), sdata, odata[0])
    return json.dumps((numpy.abs(integrate.quad(InterpolateData(isect[0], isect[1]), 0, mval)[0])/numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), 0, mval)[0])))

@require(GetData, GetNumData, InterpolateData, SampleIntersection, GetRegion)
def GraphSimilarityRegion(Samples, Primers_ID_F, Primers_ID_R, Start_Region, Stop_Region):
    svars = list(list(json.loads(j[3]) for j in sm.SampleVariables(dbs, sampleID = i, ) if (lambda x: True if "Primers_ID_F" in x.keys() and x["Primers_ID_F"] == Primers_ID_F and x["Primers_ID_R"] == Primers_ID_R else False)(json.loads(j[3])))[0] for i in Samples)
    odata = list((lambda x: (lambda y: [y[:,0], y[:,1]])(GetRegion(x[0], Start_Region, Stop_Region, x[1], True)))(GetNumData(GetData(open(svars[i]['Path'] + svars[i]['File_Prefix'] + svars[i]['Data_RC'] + svars[i]['File_Suffix']).read().splitlines()))) for i in range(len(Samples)))
    wrapprint(odata)
    maxval = max(list(max(i[0]) for i in odata))
    minval = max(list(min(i[0]) for i in odata))
    sdata = list([i[0], i[1]*(numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), minval, maxval)[0])/numpy.abs(integrate.quad(InterpolateData(i[0], i[1]), minval, maxval)[0]))] for i in odata[1:])
    isect = reduce(lambda x, y: SampleIntersection(x[0], x[1], y[0], y[1], start=minval), sdata, odata[0])
    return json.dumps((numpy.abs(integrate.quad(InterpolateData(isect[0], isect[1]), minval, maxval)[0])/numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), minval, maxval)[0])))

@require(GetData, GetNumData, InterpolateData, SampleIntersection, GetRegion)
def GraphSimilarityRegions(Samples, Primers_ID_F, Primers_ID_R, Start_Region, Stop_Region):
    svars = list(list(json.loads(j[3]) for j in sm.SampleVariables(dbs, sampleID = i, ) if (lambda x: True if "Primers_ID_F" in x.keys() and x["Primers_ID_F"] == Primers_ID_F and x["Primers_ID_R"] == Primers_ID_R else False)(json.loads(j[3]))) for i in Samples)
    odata = list(list((lambda x: (lambda y: [y[:,0], y[:,1]])(GetRegion(x[0], Start_Region, Stop_Region, x[1], True)))(GetNumData(GetData(open(j['Path'] + j['File_Prefix'] + j['Data_RC'] + j['File_Suffix']).read().splitlines()))) for j in i) for i in svars)
    maxval = max(list(max(list(max(j[0]) for j in i) for i in odata)))
    minval = min(list(min(list(min(j[0]) for j in i) for i in odata)))
    sdata = list(list([j[0], j[1] * (float(100)/numpy.abs(integrate.quad(InterpolateData(j[0], j[1]), minval, maxval)[0]))] for j in i) for i in odata)
    isect = list(list(SampleIntersection(j[0], j[1], k[0], k[1], start=minval) for j in i[0] for k in i[1]) for i in combinations(sdata, 2))[0]
    return json.dumps(list((numpy.abs(integrate.quad(InterpolateData(i[0], i[1]), minval, maxval)[0])/float(100)) for i in isect))

@require(GetData, GetNumData, InterpolateData, SampleIntersection, GetRegion, DefineRegions)
def GraphDisSimRegions(Samples, Primers_ID_F, Primers_ID_R, Start_Region, Stop_Region):
    svars = list(list(json.loads(j[3]) for j in sm.SampleVariables(dbs, sampleID = i, ) if (lambda x: True if "Primers_ID_F" in x.keys() and x["Primers_ID_F"] == Primers_ID_F and x["Primers_ID_R"] == Primers_ID_R else False)(json.loads(j[3]))) for i in Samples)
    odata = list(list((lambda x: (lambda y: [y[:,0], y[:,1]])(GetRegion(x[0], Start_Region, Stop_Region, x[1], True)))(GetNumData(GetData(open(j['Path'] + j['File_Prefix'] + j['Data_RC'] + j['File_Suffix']).read().splitlines()))) for j in i) for i in svars)
    maxval = max(list(max(list(max(j[0]) for j in i) for i in odata)))
    minval = min(list(min(list(min(j[0]) for j in i) for i in odata)))
    sdata = list(list([j[0], j[1] * (float(100)/numpy.abs(integrate.quad(InterpolateData(j[0], j[1]), minval, maxval)[0]))] for j in i) for i in odata)
    qdata =  list(list([j[0], list((j[1][q]/k[1][q] if k[1][q] >= j[1][q] else k[1][q]/j[1][q]) if k[1][q] > 0 and j[1][q] > 0 else 0 for q in range(len(j[1])))] for j in i[0] for k in i[1]) for i in combinations(sdata, 2))
    index = list([i, list([j, k] for j in range(len(sdata[i[0]])) for k in range(len(sdata[i[1]])))] for i in combinations(range(len(sdata)), 2))
    regions  = list(list(DefineRegions(sdata[index[i][0][0]][index[i][1][j][0]], sdata[index[i][0][1]][index[i][1][j][1]], qdata[i][j]) for j in range(len(index[i][1]))) for i in range(len(index)))
    return sdata, qdata, regions
    isect = list(list(SampleIntersection(j[0], j[1], k[0], k[1], start=minval) for j in i[0] for k in i[1]) for i in combinations(sdata, 2))[0]
    return json.dumps(list((numpy.abs(integrate.quad(InterpolateData(i[0], i[1]), minval, maxval)[0])/float(100)) for i in isect))

@require(GetData, GetNumData, InterpolateData, SampleIntersection)
def CorSimilarity(Samples, Primers_ID_F, Primers_ID_R):
    svars = list(list(json.loads(j[3]) for j in sm.SampleVariables(dbs, sampleID = i) if (lambda x: True if "Primers_ID_F" in x.keys() and x["Primers_ID_F"] == Primers_ID_F and x["Primers_ID_R"] == Primers_ID_R else False)(json.loads(j[3])))[0] for i in Samples)
    odata = list(GetNumData(GetData(open(svars[i]['Path'] + svars[i]['File_Prefix'] + svars[i]['Data_RC'] + svars[i]['File_Suffix']).read().splitlines())) for i in range(len(Samples)))
    mval = max(list(max(i[0]) for i in odata))
    sdata = list([i[0], i[1]*(numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), 0, mval)[0])/numpy.abs(integrate.quad(InterpolateData(i[0], i[1]), 0, mval)[0]))] for i in odata[1:])
    ccor = signal.correlate(odata[0], sdata[0])
    plt.plot(odata[0][0], odata[0][1])
    plt.plot(sdata[0][0], sdata[0][1])
    plt.plot(ccor[0], ccor[1])
    wrapprint([odata[0], sdata[0], ccor])
    print(integrate.quad(InterpolateData(ccor[0], ccor[1]), 0, mval))
    return json.dumps((numpy.abs(integrate.quad(InterpolateData(ccor[0], ccor[1]), 0, mval)[0])/numpy.abs(integrate.quad(InterpolateData(odata[0][0], odata[0][1]), 0, mval)[0])))

@require(GraphSimilarity, GraphSimilarityRegion, GraphArcDis)
def GetFunc(arow):
    global dbq, dbs, mb, sm
    mb = importlib.machinery.SourceFileLoader("mb_database", "./mb_database.py").load_module()
    sm = importlib.machinery.SourceFileLoader("sm_database", "./sm_database.py").load_module()
    dbq = mb.ConnectDatabase('./microbiome.db')
    dbs = sm.ConnectDatabase('./samples.db')
    matplotlib.use('Agg')

    if arow[2] == "GraphSimilarity":
        return [1, arow, GraphSimilarity(**dict({"Samples": json.loads(arow[1])}, **json.loads(arow[3])))]
    if arow[2] == "GraphSimilarityRegion":
        return [1, arow, GraphSimilarityRegion(**dict({"Samples": json.loads(arow[1])}, **json.loads(arow[3])))]
    if arow[2] == "GraphArcDis":
        return [1, arow, GraphArcDis(**dict({"Samples": json.loads(arow[1])}, **json.loads(arow[3])))]

def Main(cl):
    with cl[:].sync_imports():
        from scipy import signal, interpolate, stats, optimize, integrate, spatial
        from functools import reduce
        from itertools import combinations
        import json
        import numpy
        import sys
        import matplotlib
        import importlib

    global dbq, dbs
    dbq = mb.ConnectDatabase('./microbiome.db')
    dbs = sm.ConnectDatabase('./samples.db')

    proc = []
    for i in sm.Analysis(dbs, status = 'NEW'):
        for k in sorted(cl, key = lambda x: x.queue_status()['queue']):
            proc.append(k.apply_async(GetFunc, i))
            sm.Analysis(dbs, ID = i[0], status = 'PROCESSING', ins = True)
            break
    while len(proc) > 0:
        for i in proc:
            if i.ready():
                res = i.get()
                if res[0] == 0:
                    sm.Analysis(dbs, ID = res[1][0], status = 'ERROR', ins = True)
                    print(res[2])
                else:
                    sm.Analysis(dbs, ID = res[1][0], status = 'DONE', ins = True)
                    sm.Results(dbs, analysisID=res[1][0], types="JSON", value=res[2], ins=True)
                proc.remove(i)
        sleep(5)
    return proc

from scipy import signal, interpolate, stats, optimize, integrate, spatial
from functools import reduce
from itertools import combinations
import json
import numpy
import sys
import matplotlib
import importlib

dbq = mb.ConnectDatabase('./microbiome.db')
dbs = sm.ConnectDatabase('./samples.db')
