#!/usr/local/bin/python3
#Filename: spacer_finder.py
#Author: Andrew Avila (Andrew.Avila@pnnl.gov)
#Description: This script parses NCBI GBK Files and finds spacers between defined genomic elements
import os
import operator
import mb_database as mb
from Bio import SeqIO

db = mb.ConnectDatabase('./microbiome.db')

def GetSequences(fPath, fType = 'genbank'):
    print(fPath)
    return list(SeqIO.parse(open(fPath), fType))

def GetOrganism(sSeq):
    return mb.Life(db, ID = False, taxonomyID = dict(i.split(":") for i in sSeq.features[0].qualifiers['db_xref'])['taxon'], domainID = mb.Domain(db, domain = sSeq.annotations['taxonomy'][0] if len(sSeq.annotations['taxonomy']) > 0 else '', ins = True)[0][0], kingdomID = False, phylumID = mb.Phylum(db, phylum = sSeq.annotations['taxonomy'][1] if len(sSeq.annotations['taxonomy']) > 1 else '', ins = True)[0][0], classsID = mb.Classs(db, classs = sSeq.annotations['taxonomy'][2] if len(sSeq.annotations['taxonomy']) > 2 and sSeq.annotations['organism'].split(' ')[0] != sSeq.annotations['taxonomy'][2] else '', ins = True)[0][0], subclassID = False, ordersID = mb.Orders(db, order = sSeq.annotations['taxonomy'][3] if len(sSeq.annotations['taxonomy']) > 3 else '', ins = True)[0][0], suborderID = False, familyID = mb.Family(db, family = sSeq.annotations['taxonomy'][4] if len(sSeq.annotations['taxonomy']) > 4 else '', ins = True)[0][0], genusID = mb.Genus(db, genus = sSeq.annotations['taxonomy'][5] if len(sSeq.annotations['taxonomy']) > 5 else sSeq.annotations['organism'].split(' ')[0].replace('[','').replace(']',''), ins = True)[0][0], speciesID = mb.Species(db, species = ' '.join(sSeq.annotations['organism'].split(' ')[0:2]), ins = True)[0][0], strainID = mb.Strain(db, strain = sSeq.features[0].qualifiers['strain'][0] if 'strain' in sSeq.features[0].qualifiers else '', ins = True)[0][0], substrainID = mb.Substrain(db, substrain = sSeq.features[0].qualifiers['sub_strain'][0] if 'sub_strain' in sSeq.features[0].qualifiers else '', ins = True)[0][0], ins = True)

def FindSpacers(sSeq, aTerm, bTerm):
    orgID = GetOrganism(sSeq)[0][0]
    print(orgID)
    its = []
    a = list(i for i in sSeq.features if (i.type == 'rRNA' and 'product' in i.qualifiers and aTerm.upper() in i.qualifiers['product'][0].upper()))
    b = list(i for i in sSeq.features if (i.type == 'rRNA' and 'product' in i.qualifiers and bTerm.upper() in i.qualifiers['product'][0].upper()))
    if len(a) > 0 and len(b) > 0:
        a = list(i for i in a if len(list(j for j in a+b if i != j and ((j.location.start >= i.location.start and j.location.end <= i.location.end) or (j.location.end >= i.location.start and j.location.end <= i.location.end)))) == 0)
        b = list(i for i in b if len(list(j for j in a+b if i != j and ((j.location.start >= i.location.start and j.location.end <= i.location.end) or (j.location.end >= i.location.start and j.location.end <= i.location.end)))) == 0)
        c = sorted(a+b, key=operator.attrgetter('location.end'))
        for i in range(0,len(c)-1):
            if c[i].qualifiers['product'][0] != c[i+1].qualifiers['product'][0]:
                spacer = sSeq.seq[c[i].location.end + 1 : c[i+1].location.start - 1]
                if len(spacer) >= 300 and len(spacer) <= 7000:
                    mb.Sequences(db, lifeID = orgID, ncbiRefSeq = sSeq.id, description = c[i].qualifiers['product'][0], start = int(c[i].location.start), stop = int(c[i].location.end), length = len(c[i]), sequence = str(sSeq.seq[c[i].location.start : c[i].location.end]), ins = True)
                    mb.Sequences(db, lifeID = orgID, ncbiRefSeq = sSeq.id, description = c[i+1].qualifiers['product'][0], start = int(c[i+1].location.start), stop = int(c[i+1].location.end), length = len(c[i+1]), sequence = str(sSeq.seq[c[i+1].location.start : c[i+1].location.end]), ins = True)
                    mb.Sequences(db, lifeID = orgID, ncbiRefSeq = sSeq.id, description = aTerm + '_' + bTerm, start = int(c[i].location.end + 1), stop = int(c[i+1].location.start - 1), length = len(spacer), sequence = str(spacer), ins = True)

def ProcessGBKS(fPath, aTerm, bTerm):
    for root, subFolders, files in os.walk(fPath):
        list(FindSpacers(GetSequences(os.path.join(root, i))[0], aTerm, bTerm) for i in files if i.split(os.extsep, 1)[1].upper() == "GBK")
        for i in subFolders:
            ProcessGBKS(os.path.join(root, i), aTerm, bTerm)
