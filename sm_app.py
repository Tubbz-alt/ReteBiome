#!/usr/local/bin/python3
#Filename: sm_app.py
#Author: Andrew Avila (Andrew.Avila@pnnl.gov)
#Description: This script provides a web interface for use in interacting with the sample database.

from bottle import route, run, get, request, debug, static_file, default_app
from itertools import combinations
import mb_database as mb
import sm_database as sm
import json
import numpy
import time
import re

dbq = mb.ConnectDatabase('./microbiome.db')
dbs = sm.ConnectDatabase('./samples.db')

def mapRC(r, c):
    rs = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return rs[r] + str(c+1)

def demapRC(rc):
    rcs = re.findall(r"[^\W\d_]+|\d+", rc)
    if len(rcs) == 1:
        return None
    rs = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return [rs.index(rcs[0])+1, int(rcs[1])]

def genTable(arr):
    return "<table border='1'>" + "".join(list("<tr>" + "".join(list("<td>" + str(j) + "</td>" for j in i)) + "</tr>" for i in arr)) + "</table>"

def genOpt(arr):
    return "".join(list("<option value='%s'>" % str(i[0]) + str(i[1]) + "</option>" for i in arr))

@route('/list')
def listSamples():
    return genTable(list(["<a href='/sample/%s'>%s</a>" % (str(i[0]), str(i[1]))] for i in sm.Samples(dbs)))

@route('/sample/<smp>')
def sample(smp):
    with dbs:
        smpVars = dbs.execute("SELECT Variable, Types, Value FROM SAMPLE_VARIABLES LEFT JOIN VARIABLES ON Variables_ID = Variables.ID WHERE Samples_ID = ?", (smp,))
    return "Sample: %s<br><a href='/sample/%s/varadd'>Add Variable</a>%s" % (smp, smp, genTable(smpVars))

@route('/sample/<smp>/varadd')
def sampleVariableAddForm(smp):
    return "<form action='addvar'>Variable: <input name='vari'><br>Variable Type: <input name='varitype'><br>Value: <input name='varival'><br><input type='submit' value='Submit'></form><br>" + genTable([["Variable", "Variable Type"]] + list([i[1],i[2]] for i in sm.Variables(dbs))) + "<br>" + genTable([["ID", "Primer Sequence 5' -> 3'"]] + mb.Primers(dbq))

@route('/sample/<smp>/addvar')
def sampleVariableAdd(smp):
    sm.SampleVariables(dbs, sampleID = int(smp), variablesID = sm.Variables(dbs, variable = request.query.getall('vari')[0], types = request.query.getall('varitype')[0], ins = True)[0][0], value = request.query.getall('varival')[0], ins = True)
    return "Variable added to sample: %s" % smp

@route('/sample/add')
def sampleAddForm():
    return "<form action='addsmp'>Sample: <input name='smp'><br><input type='submit' value='Submit'></form>"

@route('/sample/addsmp')
def sampleAdd():
    sm.Samples(dbs, sample = request.query.getall('smp')[0], ins = True)
    return "Sample added."

@route('/find')
def findSamples():
    return "<form action='find/query'>Sample: <input name='samp'><br>Variables: <input name='vari'><br><input type='submit' value='Submit'></form>"

@route('/find/query')
def queryForm():
    samp = str(request.query.getall('samp')[0])
    vari = list(i.split('=') for i in str(request.query.getall('vari')[0]).split(',') if len(i.split('=')) > 1)
    return ""

@route('/')
def samples():
    return "<a href='/list'>List all samples</a><br><a href='/find'>Find a sample</a><br><a href='/sample/add'>Add a new sample</a><br><a href='/analysis'>List Analyses</a><br><a href='/analysis/add'>Add Analyses</a>"

@route('/results/<res>/json')
def resultsJSON(res):
    return (lambda x: {'Analysis' : x[1], 'Types' : x[2], 'Value' : x[3]})(sm.Results(dbs, ID = int(res))[0])

@route('/analysis/<id>/results/json')
def analysisResultsJSON(id):
    with dbs:
        return (lambda x: {'Samples_IDS' : json.loads(x[0]), 'Analysis_Types' : x[1], 'AVars' : json.loads(x[2]), 'Status' : x[3], 'Results_Types' : x[4], 'Value' : json.loads(x[5])})(dbs.execute('SELECT Samples_IDS, Analysis.Types, Analysis.AVars, Status, Results.Types, Value FROM Results LEFT JOIN Analysis ON Analysis_ID = Analysis.ID WHERE Analysis.ID = ?', (id,)).fetchall()[0])

@route('/analysis')
def listAnalysis():
    return genTable([["Analysis ID", "Samples", "Analysis Type", "Analysis Variables", "Analysis Status"]] + list(["<a href='/analysis/" + str(i[0]) + "/results/json'>" + str(i[0]) + "</a>"] + list(i[1:]) if i[-1] == "DONE" else i for i in sm.Analysis(dbs)))

@route('/analysis/add')
def formAnalysis():
    return "<form action='addanalysis'>Samples: <input name='samps'><br>Analysis Type: <input name='atype'><br>Analysis Variables: <input name='avars'><br><input type='submit' value='Submit'></form>"

@route('/analysis/addanalysis')
def addAnalysis():
    analysis = sm.Analysis(dbs, sampleIDS = request.query.getall('samps')[0], types = request.query.getall('atype')[0], avars = request.query.getall('avars')[0], status = 'NEW', ins = True)
    return "Analysis Added: " + analysis[0]

@route('/mb/faddsample')
def formAnalysis():
    return "<form action='addsample'>Sample: <input name='sample'><br>Sample Key: <input name='samplekey'><br>Date Sample Collected: <input name='datesample'><br>Date Sample Processed: <input name='dateprocess'><br>Part of Body Sampled: <input name='partnum'><br>Sample Time Period: <input name='sampleperiod'><br><input type='submit' value='Submit'></form>"

@route('/mb/addsample')
def mbAddSample():
    fsample = request.query.getall('sample')[0]
    fkey = request.query.getall('samplekey')[0]
    fdatesample = request.query.getall('datesample')[0]
    fdateprocess = request.query.getall('dateprocess')[0]
    fpartnum = int(request.query.getall('partnum')[0])
    fsampleperiod = int(request.query.getall('sampleperiod')[0])
    allparts = 7 if fpartnum == -1 else 2
    for i in range(1, allparts):
        fpartnum = i if allparts == 7 else fpartnum
        fpartstr = "Index Finger" if fpartnum == 1 else ("Palm" if fpartnum == 2 else ("Philtrum" if fpartnum == 3 else ("Nasal Vestibule" if fpartnum == 4 else ("Pencil" if fpartnum == 5 else ("Scissor Control" if fpartnum == 6 else "Unknown")))))
        ftube = fkey + "-" + str(fpartnum) + "-" + str(fsampleperiod)
        fdict = {"Date_Sampled" : fdatesample, "Date_Processed" : fdateprocess, "Part_Num" : fpartnum, "Part_Str" : fpartstr, "Sample_Period" : fsampleperiod, "Tube" : ftube}
    
        sid = sm.Samples(dbs, sample = fsample, ins = True)
        skey = sm.Variables(dbs, variable = "sample_key", types = "text", ins = True)
        ssample = sm.Variables(dbs, variable = "sample", types = "JSON", ins = True)
        sm.SampleVariables(dbs, sampleID = sid[0][0], variablesID = skey[0][0], value = fkey, ins = True)
        sm.SampleVariables(dbs, sampleID = sid[0][0], variablesID = ssample[0][0], value = json.dumps(fdict), ins = True)
    return "Sample added"

@route('/mb/fgenplates')
def formGenPlates():
    return "<form action='genplates'>Samples:<br><select name='samples' multiple='multiple' size=5>%s</select><br>Part Number(s): <input name='partnums'><br>Sample Point(s): <input name='sampoint'><br>Primer Pair(s): <input name='primers'><br>Number of Replicates: <input name='reps'><br>Number of Water Controls per Primer: <input name='water'><br>Number of Wells in PCR Plate: <input name='wells'><br><input type='checkbox' name='unique' value='1' checked>Generate for samples that have not yet been PCR'ed with the previous parameters<br><input type='submit' value='Submit'></form>" % (list("<option value='" + str(i[0]) +"'>" + i[1] + "</option>" for i in sm.Samples(dbs)))

@route('/mb/genplates')
def genPlates():
    getsampkey = lambda x: dbs.execute("SELECT Value FROM SAMPLE_VARIABLES LEFT JOIN VARIABLES ON VARIABLES_ID = VARIABLES.ID LEFT JOIN SAMPLES on SAMPLES_ID = SAMPLES.ID WHERE Variable = 'sample_key' AND SAMPLES.ID = ?", (x,)).fetchall()[0][0]
    fsamples = sorted(list(getsampkey(sm.Samples(dbs, ID = i.replace(' ', ''))[0][0]) for i in request.query.getall('samples')), key = lambda x: int(x[1:]))
    fpartnums = list(i.replace(' ', '') for i in request.query.getall('partnums')[0].split(','))
    fsampoint = list(i.replace(' ', '') for i in request.query.getall('sampoint')[0].split(','))
    fprimers = list(i.replace(' ', '') for i in request.query.getall('primers')[0].split(','))
    fwells = int(request.query.getall('wells')[0])
    fwater = int(request.query.getall('water')[0])
    freps = request.query.getall('reps')[0]
    funique = request.query.getall('unique')[0]
    speedvacs = []
    sampstack = []
    pstacks = []
    epstacks = []
    ssplit = lambda x: [int(x[0])] if len(x) == 1 else list(range(int(x[0]),int(x[1])+1))
    for i in fprimers:
        psplit = i.split(':')
        fsplit = sum(list(ssplit(j.split('-')) for j in psplit[0].split("'")), [])
        rsplit = sum(list(ssplit(j.split('-')) for j in psplit[1].split("'")), [])
        istack = []
        for j in fsplit:
            for k in rsplit:
                istack += [[j,k]]
        pstacks += [istack]
        epstacks += [[]]
    for i in fsamples:
        for j in fpartnums:
            for k in fsampoint:
                speedvacs += [json.loads(dbs.execute("SELECT Value FROM SAMPLE_VARIABLES LEFT JOIN VARIABLES ON VARIABLES_ID = VARIABLES.ID LEFT JOIN SAMPLES on SAMPLES_ID = SAMPLES.ID WHERE Variable = 'sample' AND SAMPLE_VARIABLES.VALUE LIKE ?", ('%' + i + '-' + j + '-' + k + '%',)).fetchall()[0][0])['Date_Processed']]
                for l in range(len(pstacks)):
                    primpair = pstacks[l].pop()
                    epstacks[l] += [primpair]
                    sampstack += [[i + '-' + j + '-' + k, primpair]]
                    if len(pstacks[l]) == 0:
                        pstacks[l] = epstacks[l][::-1]
                        epstacks[l] = []
    for i in set(speedvacs):
        for j in range(1,3):
            for l in range(len(pstacks)):
                primpair = pstacks[l].pop()
                epstacks[l] += [primpair]
                sampstack += [['S-' + '-'.join(i.split('/')) + '-' + str(j), primpair]]
                if len(pstacks[l]) == 0:
                    pstacks[l] = epstacks[l][::-1]
                    epstacks[l] = []
    for i in range(1,fwater+1):
        for l in range(len(pstacks)):
            primpair = pstacks[l].pop()
            epstacks[l] += [primpair]
            sampstack += [['W-' + '-'.join(time.strftime("%m/%d/%y").lstrip('0').split('/')) + '-' + str(i), primpair]]
            if len(pstacks[l]) == 0:
                pstacks[l] = epstacks[l][::-1]
                epstacks[l] = []
    sorttube = lambda y: [int(y[0][1:]) if y[0] != 'S' and y[0] != 'W' else 1000000 + int(y[1]), int(y[1]) if y[0] != 'S' and y[0] != 'W' else 1000000 + int(y[2]), int(y[2]) if y[0] != 'S' and y[0] != 'W' else 1000000 + int(y[3])]
    sampstack = sorted(sampstack, key = lambda x: (x[1], sorttube(x[0].split('-'))))[::-1]

    pcount = {}
    for i in sampstack:
        ppairstr = ','.join(str(j) for j in i[1])
        if ppairstr in pcount:
            pcount[ppairstr] += 1
        else:
            pcount[ppairstr] = 1
    mmix = [['Number of Reactions (Variance Added)', 'F Primer', 'R Primer', '5X OneTaq Standard Reaction Buffer', '10 mM dNTPs', '10 uM Forward Primer', '10 uM Reverse Primer', 'OneTaq DNA Polymerase', 'Template DNA', 'Nuclease Free Water', 'Amount of Master Mix'] ]
    mmix += list([str(pcount[key]) + '(' + str(int(numpy.ceil(pcount[key] * 0.05))) + ')', key.split(',')[0], key.split(',')[1], 5.0 * (pcount[key] + numpy.ceil(pcount[key] * 0.05)), 0.5 * (pcount[key] + numpy.ceil(pcount[key] * 0.05)), 0.25 * (pcount[key] + numpy.ceil(pcount[key] * 0.05)), 0.25 * (pcount[key] + numpy.ceil(pcount[key] * 0.05)), 0.125 * (pcount[key] + numpy.ceil(pcount[key] * 0.05)), 2, 17.375 * (pcount[key] + numpy.ceil(pcount[key] * 0.05)), 23] for key in sorted(pcount.keys()))

    pcrid = sm.Variables(dbs, variable = 'pcr', types = 'JSON', ins = True)[0][0]
    plates = []
    for i in range(int(numpy.ceil(len(sampstack)/fwells))):
        if fwells == 96:
            plate = numpy.zeros((8,12)).tolist()
        for j in range(len(plate)):
            for k in range(len(plate[0])):
                if len(sampstack) > 0:
                    psamp = sampstack.pop()
                    pssamp =  psamp[0].split('-')
                    if pssamp[0] != 'S' and pssamp[0] != 'W':
                        pdict = {"Plate" : time.strftime("%m/%d/%y").lstrip('0') + '-' + str(i+1), "Primers_ID_F" : psamp[1][0], "Primers_ID_R" : psamp[1][1], "Plate_RC" : mapRC(j,k), "Sample_Key" : pssamp[0], "Part_Num" : pssamp[1], "Sample_Period" : pssamp[2], "Tube": psamp[0], "Date" : time.strftime("%m/%d/%y").lstrip('0')}
                        sampid = dbs.execute("SELECT SAMPLES_ID FROM SAMPLE_VARIABLES WHERE VALUE = ?", (pssamp[0],)).fetchall()[0][0]
                        sm.SampleVariables(dbs, sampleID = sampid, variablesID = pcrid, value = json.dumps(pdict), ins = True)
                    else:
                        pdict = {"Plate" : time.strftime("%m/%d/%y").lstrip('0') + '-' + str(i+1), "Primers_ID_F" : psamp[1][0], "Primers_ID_R" : psamp[1][1], "Plate_RC" : mapRC(j,k), "Control_Type" : pssamp[0], "Control_Num" : pssamp[-1], "Tube": psamp[0], "Date" : time.strftime("%m/%d/%y").lstrip('0')}
                        sampid = dbs.execute("SELECT ID FROM SAMPLES WHERE SAMPLE = 'Human Microbiome Controls'", ()).fetchall()[0][0]
                        sm.SampleVariables(dbs, sampleID = sampid, variablesID = pcrid, value = json.dumps(pdict), ins = True)
                    plate[j][k] = psamp[0] + "<br>" + str(psamp[1][0]) + "," + str(psamp[1][1])
                else:
                    plate[j][k] = None
        plates += [plate]

    return '<br>'.join(list(time.strftime("%m/%d/%y").lstrip('0') + '-' + str(i+1) + '<br>' + genTable(plates[i]) for i in range(len(plates)))) + '<br>' + genTable(mmix) + '<br>' + genTable([['Primer ID', 'Primer Sequence', 'Primer Name']] + list(i for i in mb.Primers(dbq) if str(i[0]) in set(','.join(pcount.keys()).split(','))))

@route('/mb/fbaplates')
def fbaplates():
    return "<form action='baplates'>Samples:<br><select name='plates' multiple='multiple' size=5>%s</select><br>Number of Wells in Bioanalyzer Plate: <input name='wells'><br>Number of Channels on Pipettor: <input name='numchan'><br><input type='submit' value='Submit'></form>" % (list("<option value='" + i +"'>" + i + "</option>" for i in sorted(set(list((json.loads(j[0]))['Plate'] for j in dbs.execute("SELECT VALUE FROM SAMPLE_VARIABLES LEFT JOIN VARIABLES ON VARIABLES_ID = VARIABLES.ID WHERE VARIABLE = 'pcr'").fetchall())))))

@route('/mb/baplates')
def baplates():
    samples = sorted(list([i[0], json.loads(i[1]), i[2]] for i in dbs.execute("SELECT SAMPLES_ID, VALUE, SAMPLE_VARIABLES.ID FROM SAMPLE_VARIABLES LEFT JOIN VARIABLES ON VARIABLES_ID = VARIABLES.ID WHERE VARIABLE = 'pcr'").fetchall() if json.loads(i[1])['Plate'] in request.query.getall('plates')), key = lambda x: [x[1]['Plate'], demapRC(x[1]['Plate_RC'])])
    fwells = int(request.query.getall('wells')[0])
    fnumchan = int(request.query.getall('numchan')[0])
    baid = sm.Variables(dbs, variable = 'bioanalyzer', types = 'JSON', ins = True)[0][0]
    plates = []
    unmapped = samples
    lastlen = len(samples)
    if fwells == 96:
        plate = numpy.zeros((8,12)).tolist()
    elif fwells == 384:
        plate = numpy.zeros((16,24)).tolist()
    while (len(samples)) > 0:
        chandiv = int(len(plate)/fnumchan)
        lastlen = len(samples)
        item = samples.pop(0)
        mapped = False
        for i in range(len(plate[0])):
            for j in range(chandiv):
                if plate[j][i] == 0 and mapped == False:
                    plate[j][i] = item
                    mapped = True
                elif plate[j][i] != 0 and plate[j][i][1]['Plate'] == item[1]['Plate'] and demapRC(plate[j][i][1]['Plate_RC'])[1] == demapRC(item[1]['Plate_RC'])[1] and mapped == False:
                    plate[demapRC(item[1]['Plate_RC'])[0] * chandiv - 2 + j][i] = item
                    mapped = True
        if mapped == False:
            samples += [item]
        if lastlen == len(samples):
            plates += [plate]
            samples = sorted(samples, key = lambda x: [x[1]['Plate'], demapRC(x[1]['Plate_RC'])])
            if fwells == 96:
                plate = numpy.zeros((8,12)).tolist()
            elif fwells == 384:
                plate = numpy.zeros((16,24)).tolist()
    plates += [plate]
    for i in range(len(plates)):
        for j in range(len(plates[i])):
            for k in range(len(plates[i][j])):
                if plates[i][j][k] == 0:
                    plates[i][j][k] = None
                else:
                    sm.SampleVariables(dbs, sampleID = plates[i][j][k][0], variablesID = baid, value = json.dumps({'PCR_ID' : plates[i][j][k][2], 'PCR_Plate' : plates[i][j][k][1]['Plate'], 'Plate' : time.strftime("%m/%d/%y").lstrip('0') + '-' + str(i+1) + '-BA', 'PCR_Plate_RC' : plates[i][j][k][1]['Plate_RC'], 'Plate_RC' : mapRC(j, k), 'Date' : time.strftime("%m/%d/%y").lstrip('0'), 'File_Prefix' : time.strftime("%m%d%y").lstrip('0') + '-' + str(i+1) + '-BA_Data_', 'File_Suffix' : '.csv'}), ins = True)
                    plates[i][j][k] = plates[i][j][k][1]['Plate'] + '<br>' + plates[i][j][k][1]['Plate_RC']
    return '<br>'.join(list('Plate: ' + time.strftime("%m/%d/%y").lstrip('0') + '-' + str(i+1) + '<br>File Name Prefix: ' + time.strftime("%m%d%y").lstrip('0') + '-' + str(i+1) + '-BA_Data_' + '<br>' + genTable(plates[i]) for i in range(len(plates))))

@route('/network/build')
def buildNetwork():
    analysistype = "GraphSimilarityRegion"
    threshold = (lambda x: numpy.mean(list(x[i] for i in numpy.random.randint(0, len(x)-1, size=10000))))(list(i[0] for i in dbs.execute("SELECT avg(Value), Samples_IDS FROM Results LEFT JOIN Analysis ON Analysis_ID = Analysis.ID WHERE Analysis.Types = ? GROUP BY Samples_IDS", (analysistype,)).fetchall()))
    graph = {"nodes" : {}, "edges" : {}}
    for i in combinations((int(i) for i in request.query.samples.split(',')),2):
        with dbs:
            if json.loads(dbs.execute("SELECT Value FROM Results LEFT JOIN Analysis ON Analysis_ID = Analysis.ID WHERE Samples_IDS = ? AND Analysis.Types = ?", (json.dumps(i), analysistype)).fetchall()[0][0]) > threshold:
                graph["nodes"].update({i[0] : {}})
                graph["nodes"].update({i[1] : {}})
                if i[0] in graph["edges"]:
                    if i[1] in graph["edges"][i[0]]:
                        graph["edges"][i[0]].update({i[1] : {"length" : (graph["edges"][i[0]][i[1]]["length"] + 10*(1 - json.loads(dbs.execute("SELECT Value FROM Results LEFT JOIN Analysis ON Analysis_ID = Analysis.ID WHERE Samples_IDS = ? AND Analysis.Types = ?", (json.dumps(i), analysistype)).fetchall()[0][0]))) / 2}})
                    else:
                        graph["edges"][i[0]].update({i[1] : {"length" : 10*(1 - json.loads(dbs.execute("SELECT Value FROM Results LEFT JOIN Analysis ON Analysis_ID = Analysis.ID WHERE Samples_IDS = ? AND Analysis.Types = ?", (json.dumps(i), analysistype)).fetchall()[0][0]))}})
                else:
                    graph["edges"].update({i[0] : {i[1] : {"length" : 10*(1 - json.loads(dbs.execute("SELECT Value FROM Results LEFT JOIN Analysis ON Analysis_ID = Analysis.ID WHERE Samples_IDS = ? AND Analysis.Types = ?", (json.dumps(i), analysistype)).fetchall()[0][0]))}}})
    return graph
    
@route('/static/<filename>')
def sStatic(filename):
    return static_file(filename, root='./static/')
    
@route('/network')
def network():
    return """
<html>
<head>
</head>
<body>
<canvas id = "netview" width = "800" height = "600"></canvas>
<script src="//ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>
<script src="/static/arbor.js"></script>
<script type="text/javascript">
(function($){

  var $_GET = {};

  document.location.search.replace(/\??(?:([^=]+)=([^&]*)&?)/g, function () {
    function decode(s) {
      return decodeURIComponent(s.split("+").join(" "));
    }

    $_GET[decode(arguments[1])] = decode(arguments[2]);
  });

  SimpleRenderer = function(canvas){
    var canvas = $(canvas).get(0)
    var ctx = canvas.getContext("2d");
    var particleSystem = null

    var that = {
      init:function(system){
        particleSystem = system
        particleSystem.screenSize(canvas.width, canvas.height) 
        particleSystem.screenPadding(80)
      },
      
      redraw:function(){
        ctx.clearRect(0,0, canvas.width, canvas.height)
        
        particleSystem.eachEdge(function(edge, pt1, pt2){
          ctx.strokeStyle = "rgba(0,0,0, .333)"
          ctx.lineWidth = 1 + 10/edge.length
          ctx.beginPath()
          ctx.moveTo(pt1.x, pt1.y)
          ctx.lineTo(pt2.x, pt2.y)
          ctx.stroke()
        })

        particleSystem.eachNode(function(node, pt){
          var w = 10
          ctx.fillStyle = "red"
          ctx.font = "bold 14px Arial"
          ctx.fillText(node.name, pt.x-w/2, pt.y-w/2)
        })    			
      }
    }
    return that
  }    

  $(document).ready(function(){
    var sys = arbor.ParticleSystem(1000, 800, 0.5)
    sys.renderer = SimpleRenderer("#netview")
    var data = $.getJSON("/network/build?samples=" + $_GET["samples"],function(data){sys.graft({nodes:data.nodes, edges:data.edges})}) 
  })
})(this.jQuery)
</script>
</body>
</html>
"""

if __name__ == '__main__':
    debug(True)
    run(reloader=True, port=8080)
else:
    os.chdir(os.path.dirname(__file__))
    application = default_app()
