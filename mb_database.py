#!/usr/local/bin/python3
#Filename: mb_database.py
#Author: Andrew Avila (Andrew.Avila@pnnl.gov)
#Description: Database layer for microbiome project
import sqlite3 as sql

def ConnectDatabase(dbPath):
    return sql.connect(dbPath)

def CreateDatabase(db):
    with db:        
        db.execute('CREATE TABLE LIFE (ID INTEGER PRIMARY KEY, NCBI_Taxonomy_ID INTEGER, Domain_ID INTEGER, Kingdom_ID INTEGER, Phylum_ID INTEGER, Classs_ID INTEGER, Subclass_ID INTEGER, Orders_ID INTEGER, Suborder_ID INTEGER, Family_ID INTEGER, Genus_ID INTEGER, Species_ID INTEGER, Strain_ID INTEGER, Substrain_ID INTEGER)')
        db.execute('CREATE TABLE SEQUENCES (ID INTEGER PRIMARY KEY, Life_ID INTEGER, NCBI_RefSeq TEXT, Description TEXT, Start INTEGER, Stop INTEGER, Length INTEGER, Sequence TEXT)')
        db.execute('CREATE TABLE DOMAIN (ID INTEGER PRIMARY KEY, _Domain TEXT)')
        db.execute('CREATE TABLE KINGDOM (ID INTEGER PRIMARY KEY, _Kingdom TEXT)')
        db.execute('CREATE TABLE PHYLUM (ID INTEGER PRIMARY KEY, _Phylum TEXT)')
        db.execute('CREATE TABLE CLASSS (ID INTEGER PRIMARY KEY, _Classs TEXT)')
        db.execute('CREATE TABLE SUBCLASS (ID INTEGER PRIMARY KEY, _Subclass TEXT)')
        db.execute('CREATE TABLE ORDERS (ID INTEGER PRIMARY KEY, _Order TEXT)')
        db.execute('CREATE TABLE SUBORDER (ID INTEGER PRIMARY KEY, _Suborder TEXT)')
        db.execute('CREATE TABLE FAMILY (ID INTEGER PRIMARY KEY, _Family TEXT)')
        db.execute('CREATE TABLE GENUS (ID INTEGER PRIMARY KEY, _Genus TEXT)')
        db.execute('CREATE TABLE SPECIES (ID INTEGER PRIMARY KEY, _Species TEXT)')
        db.execute('CREATE TABLE STRAIN (ID INTEGER PRIMARY KEY, _Strain TEXT)')
        db.execute('CREATE TABLE SUBSTRAIN (ID INTEGER PRIMARY KEY, _Substrain TEXT)')
        db.execute('CREATE TABLE PRIMERS (ID INTEGER PRIMARY KEY, PSequence TEXT, PName TEXT)')
        db.execute('CREATE TABLE PRIMERS_LIFE (ID INTEGER PRIMARY KEY, Primers_ID INTEGER, Life_ID INTEGER)')
        db.execute('CREATE TABLE PRIMERS_SEQUENCES (ID INTEGER PRIMARY KEY, Primers_ID INTEGER, Sequences_ID INTEGER, StartLoc INTEGER, Mismatches INTEGER)')
        db.execute('CREATE TABLE PRIMERS_REGION (ID INTEGER PRIMARY KEY, Primers_Sequences_ID_A INTEGER, Primers_Sequences_ID_B INTEGER, Region INTEGER)')

def Life(db, ID = False, taxonomyID = False, domainID = False, kingdomID = False, phylumID = False, classsID = False, subclassID = False, ordersID = False, suborderID = False, familyID = False, genusID = False, speciesID = False, strainID = False, substrainID = False, ins = False):
    with db:
        if ID == False and taxonomyID == False and domainID == False and kingdomID == False and phylumID == False and classsID == False and subclassID == False and ordersID == False and suborderID == False and familyID == False and genusID == False and speciesID == False and strainID == False and substrainID == False:
            return db.execute('SELECT ID, NCBI_Taxonomy_ID, Domain_ID, Kingdom_ID, Phylum_ID, Classs_ID, Subclass_ID, Orders_ID, Suborder_ID, Family_ID, Genus_ID, Species_ID, Strain_ID, Substrain_ID FROM LIFE').fetchall()
        elif ID == False and taxonomyID != False and domainID != False and kingdomID != False and phylumID != False and classsID != False and subclassID != False and ordersID != False and suborderID != False and familyID != False and genusID != False and speciesID != False and strainID != False and substrainID != False:
            return (lambda x: [(db.execute('INSERT INTO LIFE (NCBI_Taxonomy_ID, Domain_ID, Kingdom_ID, Phylum_ID, Classs_ID, Subclass_ID, Orders_ID, Suborder_ID, Family_ID, Genus_ID, Species_ID, Strain_ID, Substrain_ID) values (?,?,?,?,?,?,?,?,?,?,?,?,?)', (taxonomyID, domainID, kingdomID, phylumID, classsID, subclassID, ordersID, suborderID, familyID, genusID, speciesID, strainID, substrainID)).lastrowid, taxonomyID, domainID, kingdomID, phylumID, classsID, subclassID, ordersID, suborderID, familyID, genusID, speciesID, strainID, substrainID)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, NCBI_Taxonomy_ID, Domain_ID, Kingdom_ID, Phylum_ID, Classs_ID, Subclass_ID, Orders_ID, Suborder_ID, Family_ID, Genus_ID, Species_ID, Strain_ID, Substrain_ID FROM LIFE WHERE NCBI_Taxonomy_ID= ? and Domain_ID= ? and Kingdom_ID= ? and Phylum_ID= ? and Classs_ID= ? and Subclass_ID= ? and Orders_ID= ? and Suborder_ID= ? and Family_ID= ? and Genus_ID= ? and Species_ID= ? and Strain_ID= ? and Substrain_ID = ?', (taxonomyID, domainID, kingdomID, phylumID, classsID, subclassID, ordersID, suborderID, familyID, genusID, speciesID, strainID, substrainID)).fetchall())
        elif ID == False and taxonomyID != False and domainID != False and kingdomID == False and phylumID != False and classsID != False and subclassID == False and ordersID != False and suborderID == False and familyID != False and genusID != False and speciesID != False and strainID != False and substrainID != False:
            return (lambda x: [(db.execute('INSERT INTO LIFE (NCBI_Taxonomy_ID, Domain_ID, Phylum_ID, Classs_ID, Orders_ID, Family_ID, Genus_ID, Species_ID, Strain_ID, Substrain_ID) values (?,?,?,?,?,?,?,?,?,?)', (taxonomyID, domainID, phylumID, classsID, ordersID, familyID, genusID, speciesID, strainID, substrainID)).lastrowid, int(taxonomyID), domainID, None, phylumID, classsID, None, ordersID, None, familyID, genusID, speciesID, strainID, substrainID)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, NCBI_Taxonomy_ID, Domain_ID, Kingdom_ID, Phylum_ID, Classs_ID, Subclass_ID, Orders_ID, Suborder_ID, Family_ID, Genus_ID, Species_ID, Strain_ID, Substrain_ID FROM LIFE WHERE NCBI_Taxonomy_ID= ? and Domain_ID= ? and Phylum_ID= ? and Classs_ID= ? and Orders_ID= ? and Family_ID= ? and Genus_ID= ? and Species_ID= ? and Strain_ID= ? and Substrain_ID = ?', (taxonomyID, domainID, phylumID, classsID, ordersID, familyID, genusID, speciesID, strainID, substrainID)).fetchall())
        if ID != False and taxonomyID == False and domainID == False and kingdomID == False and phylumID == False and classsID == False and subclassID == False and ordersID == False and suborderID == False and familyID == False and genusID == False and speciesID == False and strainID == False and substrainID == False:
            return db.execute('SELECT ID, NCBI_Taxonomy_ID, Domain_ID, Kingdom_ID, Phylum_ID, Classs_ID, Subclass_ID, Orders_ID, Suborder_ID, Family_ID, Genus_ID, Species_ID, Strain_ID, Substrain_ID FROM LIFE WHERE ID = ?', (ID,)).fetchall()
        
def Sequences(db, ID = False, lifeID = False, ncbiRefSeq = False, description = False, start = False, stop = False, length = False, sequence = False, ins = False):
    with db:
        if ID == False and lifeID == False and ncbiRefSeq == False and description == False and start == False and stop == False and length == False and sequence == False:
            return db.execute('SELECT ID, Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence FROM SEQUENCES').fetchall()
        elif ID == False and lifeID != False and ncbiRefSeq != False and description != False and start != False and stop != False and length != False and sequence != False:
            return (lambda x: [(db.execute('INSERT INTO SEQUENCES (Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence) values (?,?,?,?,?,?,?)', (lifeID, ncbiRefSeq, description, start, stop, length, sequence)).lastrowid, lifeID, ncbiRefSeq, description, start, stop, length, sequence)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence FROM SEQUENCES WHERE Life_ID = ? and NCBI_RefSeq = ? and Description = ? and Start = ? and Stop = ?', (lifeID, ncbiRefSeq, description, start, stop)).fetchall())
        elif ID == False and lifeID == False and ncbiRefSeq == False and description == False and start == False and stop == False and length == False and sequence == False:
            return db.execute('SELECT ID, Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence FROM SEQUENCES WHERE Description = ?', (description,)).fetchall()
        elif ID == False and lifeID != False and ncbiRefSeq == False and description == False and start == False and stop == False and length == False and sequence == False:
            return db.execute('SELECT ID, Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence FROM SEQUENCES WHERE Life_ID = ?', (lifeID,)).fetchall()
        elif ID == False and lifeID != False and ncbiRefSeq == False and description == False and start != False and stop == False and length == False and sequence == False:
            return db.execute('SELECT ID, Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence FROM SEQUENCES WHERE Life_ID = ? AND Start = ?', (lifeID,start)).fetchall()
        elif ID == False and lifeID != False and ncbiRefSeq == False and description == False and start == False and stop != False and length == False and sequence == False:
            return db.execute('SELECT ID, Life_ID, NCBI_RefSeq, Description, Start, Stop, Length, Sequence FROM SEQUENCES WHERE Life_ID = ? AND Stop = ?', (lifeID,stop)).fetchall()
                
def Substrain(db, ID = False, substrain = False, ins = False):
    with db:
        if ID == False and substrain == False:
            return db.execute('SELECT ID, _Substrain FROM SUBSTRAIN').fetchall()
        elif ID == False and substrain != False:
            return (lambda x: [(db.execute('INSERT INTO SUBSTRAIN (_Substrain) values (?)', (substrain,)).lastrowid, substrain)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Substrain FROM SUBSTRAIN WHERE _Substrain = ?', (substrain,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Substrain FROM SUBSTRAIN WHERE ID = ?', (ID,)).fetchall()

def Strain(db, ID = False, strain = False, ins = False):
    with db:
        if ID == False and strain == False:
            return db.execute('SELECT ID, _Strain FROM STRAIN').fetchall()
        elif ID == False and strain != False:
            return (lambda x: [(db.execute('INSERT INTO STRAIN (_Strain) values (?)', (strain,)).lastrowid, strain)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Strain FROM STRAIN WHERE _Strain = ?', (strain,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Strain FROM STRAIN WHERE ID = ?', (ID,)).fetchall()

def Species(db, ID = False, species = False, ins = False):
    with db:
        if ID == False and species == False:
            return db.execute('SELECT ID, _Species FROM SPECIES').fetchall()
        elif ID == False and species != False:
            return (lambda x: [(db.execute('INSERT INTO SPECIES (_Species) values (?)', (species,)).lastrowid, species)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Species FROM SPECIES WHERE _Species = ?', (species,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Species FROM SPECIES WHERE ID = ?', (ID,)).fetchall()

def Genus(db, ID = False, genus = False, ins = False):
    with db:
        if ID == False and genus == False:
            return db.execute('SELECT ID, _Genus FROM GENUS').fetchall()
        elif ID == False and genus != False:
            return (lambda x: [(db.execute('INSERT INTO GENUS (_Genus) values (?)', (genus,)).lastrowid, genus)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Genus FROM GENUS WHERE _Genus = ?', (genus,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Genus FROM GENUS WHERE ID = ?', (ID,)).fetchall()

def Family(db, ID = False, family = False, ins = False):
    with db:
        if ID == False and family == False:
            return db.execute('SELECT ID, _Family FROM FAMILY').fetchall()
        elif ID == False and family != False:
            return (lambda x: [(db.execute('INSERT INTO FAMILY (_Family) values (?)', (family,)).lastrowid, family)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Family FROM FAMILY WHERE _Family = ?', (family,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Family FROM FAMILY WHERE ID = ?', (ID,)).fetchall()

def Suborder(db, ID = False, suborder = False, ins = False):
    with db:
        if ID == False and suborder == False:
            return db.execute('SELECT ID, _Suborder FROM SUBORDER').fetchall()
        elif ID == False and suborder != False:
            return (lambda x: [(db.execute('INSERT INTO SUBORDER (_Suborder) values (?)', (suborder,)).lastrowid, suborder)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Suborder FROM SUBORDER WHERE _Suborder = ?', (suborder,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Suborder FROM SUBORDER WHERE ID = ?', (ID,)).fetchall()

def Orders(db, ID = False, order = False, ins = False):
    with db:
        if ID == False and order == False:
            return db.execute('SELECT ID, _Order FROM ORDERS').fetchall()
        elif ID == False and order != False:
            return (lambda x: [(db.execute('INSERT INTO ORDERS (_Order) values (?)', (order,)).lastrowid, order)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Order FROM ORDERS WHERE _Order = ?', (order,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Order FROM ORDERS WHERE ID = ?', (ID,)).fetchall()

def Subclass(db, ID = False, subclass = False, ins = False):
    with db:
        if ID == False and subclass == False:
            return db.execute('SELECT ID, _Subclass FROM SUBCLASS').fetchall()
        elif ID == False and subclass != False:
            return (lambda x: [(db.execute('INSERT INTO SUBCLASS (_Subclass) values (?)', (subclass,)).lastrowid, subclass)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Subclass FROM SUBCLASS WHERE _Subclass = ?', (subclass,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Subclass FROM SUBCLASS WHERE ID = ?', (ID,)).fetchall()

def Classs(db, ID = False, classs = False, ins = False):
    with db:
        if ID == False and classs == False:
            return db.execute('SELECT ID, _Classs FROM CLASSS').fetchall()
        elif ID == False and classs != False:
            return (lambda x: [(db.execute('INSERT INTO CLASSS (_Classs) values (?)', (classs,)).lastrowid, classs)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Classs FROM CLASSS WHERE _Classs = ?', (classs,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Classs FROM CLASSS WHERE ID = ?', (ID,)).fetchall()

def Phylum(db, ID = False, phylum = False, ins = False):
    with db:
        if ID == False and phylum == False:
            return db.execute('SELECT ID, _Phylum FROM PHYLUM').fetchall()
        elif ID == False and phylum != False:
            return (lambda x: [(db.execute('INSERT INTO PHYLUM (_Phylum) values (?)', (phylum,)).lastrowid, phylum)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Phylum FROM PHYLUM WHERE _Phylum = ?', (phylum,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Phylum FROM PHYLUM WHERE ID = ?', (ID,)).fetchall()

def Kingdom(db, ID = False, kingdom = False, ins = False):
    with db:
        if ID == False and kingdom == False:
            return db.execute('SELECT ID, _Kingdom FROM KINGDOM').fetchall()
        elif ID == False and kingdom != False:
            return (lambda x: [(db.execute('INSERT INTO KINGDOM (_Kingdom) values (?)', (kingdom,)).lastrowid, kingdom)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Kingdom FROM KINGDOM WHERE _Kingdom = ?', (kingdom,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Kingdom FROM KINGDOM WHERE ID = ?', (ID,)).fetchall()

def Domain(db, ID = False, domain = False, ins = False):
    with db:
        if ID == False and domain == False:
            return db.execute('SELECT ID, _Domain FROM DOMAIN').fetchall()
        elif ID == False and domain != False:
            return (lambda x: [(db.execute('INSERT INTO DOMAIN (_Domain) values (?)', (domain,)).lastrowid, domain)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, _Domain FROM DOMAIN WHERE _Domain = ?', (domain,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, _Domain FROM DOMAIN WHERE ID = ?', (ID,)).fetchall()

def Primers(db, ID = False, psequence = False, pname = False, ins = False):
    with db:
        if ID == False and psequence == False and pname == False:
            return db.execute('SELECT ID, PSequence, PName FROM PRIMERS').fetchall()
        elif ID == False and psequence != False and pname != False:
            return (lambda x: [(db.execute('INSERT INTO PRIMERS (PSequence, PName) values (?,?)', (psequence,pname)).lastrowid, psequence)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, PSequence, PName FROM PRIMERS WHERE PSequence = ?', (psequence,)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, PSequence, PName FROM PRIMERS WHERE ID = ?', (ID,)).fetchall()

def PrimersLife(db, ID = False, primersID = False, lifeID = False, ins = False):
    with db:
        if ID == False and primersID == False and lifeID == False:
            return db.execute('SELECT ID, Primers_ID, Life_ID FROM PRIMERS_LIFE').fetchall()
        elif ID == False and primersID != False and lifeID != False:
            return (lambda x: [(db.execute('INSERT INTO PRIMERS_LIFE (Primers_ID, Life_ID) values (?,?)', (primersID,lifeID)).lastrowid, primersID, lifeID)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Primers_ID, Life_ID FROM PRIMERS_LIFE WHERE Primers_ID = ? and Life_ID = ?', (primersID,lifeID)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, Primers_ID, Life_ID FROM PRIMERS_LIFE WHERE ID = ?', (ID,)).fetchall()

def PrimersSequences(db, ID = False, primersID = False, sequencesID = False, startloc = False, mismatches = False, ins = False):
    with db:
        if ID == False and primersID == False and sequencesID == False and startloc == False:
            return db.execute('SELECT ID, Primers_ID, Sequences_Id, StartLoc, Mismatches FROM PRIMERS_SEQUENCES').fetchall()
        elif ID == False and primersID != False and sequencesID != False and startloc != False and mismatches != False:
            return (lambda x: [(db.execute('INSERT INTO PRIMERS_SEQUENCES (Primers_ID, Sequences_Id, StartLoc, Mismatches) values (?,?,?,?)', (primersID,sequencesID,startloc,mismatches)).lastrowid, primersID, sequencesID, startloc,mismatches)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Primers_ID, Sequences_Id, StartLoc, Mismatches FROM PRIMERS_SEQUENCES WHERE Primers_ID = ? and Sequences_Id = ? and StartLoc = ? and Mismatches = ?', (primersID,sequencesID,startloc,mismatches)).fetchall())
        elif ID == False and primersID != False and sequencesID == False and startloc == False:
            return db.execute('SELECT ID, Primers_ID, Sequences_Id, StartLoc, Mismatches FROM PRIMERS_SEQUENCES WHERE Primers_ID = ?', (primersID,)).fetchall()
        elif ID != False:
            return db.execute('SELECT ID, Primers_ID, Sequences_Id, StartLoc, Mismatches FROM PRIMERS_SEQUENCES WHERE ID = ?', (ID,)).fetchall()

def PrimersRegion(db, ID = False, primerA = False, primerB = False, region = False, ins = False):
    with db:
        if ID == False and primerA == False and primerB == False and region == False:
            return db.execute('SELECT ID, Primers_Sequences_ID_A, Primers_Sequences_ID_B, Region FROM PRIMERS_REGION').fetchall()
        elif ID == False and primerA != False and primerB != False and region != False:
            return (lambda x: [(db.execute('INSERT INTO PRIMERS_REGION (Primers_Sequences_ID_A, Primers_Sequences_ID_B, Region) values (?,?,?)', (primerA,primerB,region)).lastrowid, primerA, primerB, region)] if len(x) == 0 and ins == True else x)(db.execute('SELECT ID, Primers_Sequences_ID_A, Primers_Sequences_ID_B, Region FROM PRIMERS_REGION WHERE Primers_Sequences_ID_A = ? and Primers_Sequences_ID_B = ? and Region = ?', (primerA,primerB,region)).fetchall())
        elif ID != False:
            return db.execute('SELECT ID, Primers_Sequences_ID_A, Primers_Sequences_ID_B, Region FROM PRIMERS_REGION WHERE ID = ?', (ID,)).fetchall()
