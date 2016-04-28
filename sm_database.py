#!/usr/local/bin/python3
#Filename: sm_database.py
#Author: Andrew Avila (Andrew.Avila@pnnl.gov)
#Description: Database layer for samples and analysis results.
import sqlite3 as sql

def ConnectDatabase(dbPath):
    return sql.connect(dbPath)

def CreateDatabase(db):
    with db:        
        db.execute('CREATE TABLE SAMPLES (ID INTEGER PRIMARY KEY, Sample TEXT)')
        db.execute('CREATE TABLE VARIABLES (ID INTEGER PRIMARY KEY, Variable TEXT, Types TEXT)')
        db.execute('CREATE TABLE SAMPLE_VARIABLES (ID INTEGER PRIMARY KEY, Samples_ID INTEGER, Variables_ID INTEGER, Value TEXT)')
        db.execute('CREATE TABLE ANALYSIS (ID INTEGER PRIMARY KEY, Samples_IDS TEXT, Types TEXT, AVars TEXT, Status TEXT)')
        db.execute('CREATE TABLE RESULTS (ID INTEGER PRIMARY KEY, Analysis_ID INTEGER, Types TEXT, Value TEXT)')
        
def Samples(db, ID = False, sample = False, ins = False):
    with db:
        if ID == False and sample == False:
            return db.execute('SELECT ID, Sample FROM SAMPLES').fetchall()
        elif ID == False and sample != False:
            return (lambda x: [(db.execute('INSERT INTO SAMPLES (Sample) values (?)', (sample,)).lastrowid, sample)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Sample FROM SAMPLES WHERE Sample = ?', (sample,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, Sample FROM SAMPLES WHERE ID = ?', (ID,)).fetchall()

def Variables(db, ID = False, variable = False, types = False, ins = False):
    with db:
        if ID == False and variable == False and types == False:
            return db.execute('SELECT ID, Variable, Types FROM VARIABLES').fetchall()
        elif ID == False and variable != False and types != False:
            return (lambda x: [(db.execute('INSERT INTO VARIABLES (Variable, Types) values (?,?)', (variable, types)).lastrowid, variable, types)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Variable, Types FROM VARIABLES WHERE Variable = ? and Types = ?', (variable, types)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, Variable, Types FROM VARIABLES WHERE ID = ?', (ID,)).fetchall()

def SampleVariables(db, ID = False, sampleID = False, variablesID = False, value = False, ins = False):
    with db:
        if ID == False and sampleID == False and variablesID == False and value == False:
            return db.execute('SELECT ID, Samples_ID, Variables_ID, Value FROM SAMPLE_VARIABLES').fetchall()
        elif ID == False and sampleID != False and variablesID != False and value != False:
            return (lambda x: [(db.execute('INSERT INTO SAMPLE_VARIABLES (Samples_ID, Variables_ID, Value) values (?,?,?)', (sampleID, variablesID, value)).lastrowid, sampleID, variablesID, value)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Samples_ID, Variables_ID, Value FROM SAMPLE_VARIABLES WHERE Samples_ID = ? and Variables_ID = ? and Value = ?', (sampleID, variablesID, value)).fetchall())
        elif ID == False and sampleID != False and variablesID == False and value == False:
            return db.execute('SELECT ID, Samples_ID, Variables_ID, Value FROM SAMPLE_VARIABLES WHERE Samples_ID = ?', (sampleID,)).fetchall()
        elif ID != False:
            return db.execute('SELECT ID, Samples_ID, Variables_ID, Value FROM SAMPLE_VARIABLES WHERE ID = ?', (ID,)).fetchall()

def Analysis(db, ID = False, sampleIDS = False, types = False, avars = False, status = False, ins = False):
    with db:
        if ID == False and sampleIDS == False and types == False and avars == False and status == False:
            return db.execute('SELECT ID, Samples_IDS, Types, AVars, Status FROM ANALYSIS').fetchall()
        elif ID == False and sampleIDS != False and types != False and avars != False and status != False:
            return (lambda x: [(db.execute('INSERT INTO ANALYSIS (Samples_IDS, Types, AVars, Status) values (?,?,?,?)', (sampleIDS, types, avars, status)).lastrowid, sampleIDS, types, avars, status)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Samples_IDS, Types, AVars, Status FROM ANALYSIS WHERE Samples_IDS = ? and Types = ? and AVars = ? and Status = ?', (sampleIDS, types, avars, status)).fetchall())
        elif ID == False and sampleIDS == False and types == False and avars == False and status != False:
            return db.execute('SELECT ID, Samples_IDS, Types, AVars, Status FROM ANALYSIS WHERE Status = ?', (status,)).fetchall()
        elif ID != False and sampleIDS == False and types == False and avars == False and status != False:
            if ins == True:
                db.execute('UPDATE ANALYSIS SET STATUS = ? WHERE ID = ?', (status, ID)) 
            return db.execute('SELECT ID, Samples_IDS, Types, AVars, Status FROM ANALYSIS WHERE ID = ? AND STATUS = ?', (ID,status)).fetchall()
        elif ID != False:
            return db.execute('SELECT ID, Samples_IDS, Types, AVars, Status FROM ANALYSIS WHERE ID = ?', (ID,)).fetchall()

def Results(db, ID = False, analysisID = False, types = False, value = False, ins = False):
    with db:
        if ID == False and analysisID == False and types == False and value == False:
            return db.execute('SELECT ID, Analysis_ID, Types, Value FROM RESULTS').fetchall()
        elif ID == False and analysisID != False and types != False and value != False:
            return (lambda x: [(db.execute('INSERT INTO RESULTS (Analysis_ID, Types, Value) values (?,?,?)', (analysisID, types, value)).lastrowid, analysisID, types, value)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Analysis_ID, Types, Value FROM RESULTS WHERE Analysis_ID = ? and Types = ? and Value = ?', (analysisID, types, value)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, Analysis_ID, Types, Value FROM RESULTS WHERE ID = ?', (ID,)).fetchall()
