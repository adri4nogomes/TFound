from PositionMatrix import PositionMatrix as pm
import sqlite3
import os
import pandas as pd

def getTFs(tfs=None, EC=None, dimer=False,Database=None):
    return getAll(target=tfs, EC=EC, dimer=dimer, Database=Database)[:len(tfs)] if tfs!=None else None

def _adjust(l):
    if(l!=None):
        l = l if isinstance(l, list) else [l]
        aux = []
        for r in l:
            try:
                aux.append(int(r))
            except:
                aux.append(r.upper())
        return aux
    return []

def getAll(target=None, remove=None, EC="All", dimer=False, removeParalogs=False, all=True, Database=None):
    target = _adjust(target)
    remove = _adjust(remove)

    #Criar dict de targets
    if(target!=[]):
        targets = {i:[] for i in target}

    if(removeParalogs):
        paralogs = getParalogs(target, type="MotifID", dimer=dimer, EC=EC, all=False)
        if(paralogs!=None):
            remove.extend(paralogs)

    db = "TF.db"
    if(os.path.isfile(db)):
        conn = sqlite3.connect(db)
        c = conn.cursor()

        query = '''SELECT DISTINCT f.absPath, TFs.StandardName, fTF.SystematicName, fTF.MotifID, fTF.SubMotif, fTF.ExpertConfidence, fTF.Dimer, fTF.fileName FROM Files f, FileTF fTF, TFs
            WHERE f.Name LIKE fTF.fileName AND (fTF.SystematicName = TFs.SystematicName OR fTF.SystematicName = TFs.StandardName)'''
        if(EC in ['High','Medium','Low']):
            query += " AND ExpertConfidence LIKE '{}'".format(EC)
        if(Database!=None):
            query += " AND f.Database LIKE '{}'".format(Database)

        result = []
        for reg in c.execute(query):
            if(reg[1] not in remove and reg[2] not in remove and reg[3] not in remove):
                absPath = reg[0]
                StandardName = reg[1] if not reg[6] else reg[7].split("_")[0]
                SystematicName = reg[2] if not reg[6] else reg[7].split("_")[0]
                MotifID = int(reg[3]) if reg[3]!=None else None
                SubMotif = int(reg[4]) if reg[4]!=None else None
                ExpertConfidence = reg[5]

                aux = pm(absPath, StandardName, SystematicName, MotifID, SubMotif, EC=ExpertConfidence)
                if((reg[1].upper() if dimer else StandardName.upper()) in target):
                    targets[reg[1].upper()].append(aux)
                elif((reg[2].upper() if dimer else SystematicName.upper()) in target):
                    targets[reg[2].upper()].append(aux)
                elif(MotifID in target):
                    targets[reg[3]].append(aux)
                else:
                    result.append(aux)
        conn.close()

        if 'targets' in locals():
            for value in reversed(list(targets.values())):
                if((not all and value!=[]) or all):
                    result.insert(0,value if len(value)>1 else (value[0] if len(value)==1 else None))

        return result
    else:
        return None

def getParalogs(tfs, type="MotifID", dimer=False, EC="All", all=True, Database=None):
    if(tfs==None or tfs==[]):
        return []
    target = _adjust(tfs)

    Names = []
    MotifIDs = []
    for tf in target:
        try:
            MotifIDs.append(int(tf))
        except:
            Names.append(tf)

    #Criar dict de targets
    if(target!=[]):
        targets = {i:[] for i in target}

    db = "TF.db"
    if(os.path.isfile(db)):
        conn = sqlite3.connect(db)
        c = conn.cursor()
        query = "SELECT DISTINCT ta.StandardName, ta.SystematicName, fa.Dimer, fa.fileName, fa.MotifID, tb.StandardName, tb.SystematicName, fb.Dimer, fb.fileName, fb.MotifID FROM Paralogs p, TFs ta, TFs tb, FileTF fa, FileTF fb, Files f "
        query += "WHERE ta.SystematicName=p.ParalogA AND tb.SystematicName=p.ParalogB AND fa.SystematicName=p.ParalogA AND fb.SystematicName=p.ParalogB AND f.Name LIKE fa.fileName AND ("
        query += 'ta.SystematicName IN ({}) OR tb.SystematicName IN ({}) OR ta.StandardName IN ({}) OR tb.StandardName IN ({})' if len(Names)>0 else ''
        query += ' OR ' if len(MotifIDs)>0 and len(Names)>0 else ''
        query += 'fa.MotifID IN({}) OR fb.MotifID IN({})' if len(MotifIDs)>0 else ''
        query +=')'
        if(EC in ['High','Medium','Low']):
            query += " AND fa.ExpertConfidence LIKE '{}' AND fb.ExpertConfidence LIKE '{}'".format(EC,EC)
        if(Database!=None):
            query += " AND f.Database LIKE '{}'".format(Database)

        a = tuple(Names)+tuple(Names)+tuple(Names)+tuple(Names)
        b = tuple(MotifIDs)+tuple(MotifIDs)
        m = ','.join(["?"]*len(MotifIDs))
        n = ','.join(["?"]*len(Names))

        if(len(MotifIDs)>0 and len(Names)>0):
            query = query.format(n, n, n, n, m, m)
            t = a+b
        elif(len(MotifIDs)>0):
            query = query.format(m,m)
            t = b
        else:
            query = query.format(n,n,n,n)
            t = a

        result = []

        for reg in c.execute(query, t):
            StandardNameA = reg[0] if not reg[2] else reg[3].split("_")[0]
            SystematicNameA = reg[1] if not reg[2] else reg[3].split("_")[0]
            MotifIDA = int(reg[4])

            StandardNameB = reg[5] if not reg[7] else reg[8].split("_")[0]
            SystematicNameB = reg[6] if not reg[7] else reg[8].split("_")[0]
            MotifIDB = int(reg[9])

            auxA = StandardNameA if(type=="StandardName") else (SystematicNameA if(type=="SystematicName") else MotifIDA if(type=="MotifID") else "{} ({})".format(StandardNameA,MotifIDA))
            auxB = StandardNameB if(type=="StandardName") else (SystematicNameB if(type=="SystematicName") else MotifIDB if(type=="MotifID") else "{} ({})".format(StandardNameB,MotifIDB))

            if((reg[0] if dimer else StandardNameA) in target):
                targets[reg[0]].append(auxB)
            elif((reg[1] if dimer else SystematicNameA) in target):
                targets[reg[1]].append(auxB)
            elif(reg[4] in target):
                targets[reg[4]].append(auxB)

            if((reg[5] if dimer else StandardNameB) in target):
                targets[reg[5]].append(auxA)
            elif((reg[6] if dimer else SystematicNameB) in target):
                targets[reg[6]].append(auxA)
            elif(reg[9] in target):
                targets[reg[9]].append(auxA)

        conn.close()

        for value in reversed(list(targets.values())):
            if((not all and value!=None) or all):
                result.insert(0,value if len(value)>1 else (value[0] if len(value)==1 else None))

        return result

def getDBs():
    db = "TF.db"
    if(os.path.isfile(db)):
        conn = sqlite3.connect(db)
        query = "SELECT DISTINCT ExpertConfidence, Database FROM FileTF ftf, Files f WHERE ftf.fileName==f.Name ORDER BY Database DESC"
        ecs = pd.read_sql_query(query,conn)
        conn.close()
        return ecs
    return None

def getRegulators():
    db = "TF.db"
    if(os.path.isfile(db)):
        conn = sqlite3.connect(db)
        query = "SELECT DISTINCT RegulatorStd, RegulatorSys FROM Regulation"
        regulators = pd.read_sql_query(query,conn)
        conn.close()
        return sorted([regulators.iloc[i,0]+"/"+regulators.iloc[i,1] for i in range(len(regulators))])
    return None

def getTargets(regulator):
    db = "TF.db"
    if(os.path.isfile(db)):
        conn = sqlite3.connect(db)
        query = "SELECT TargetSys, TargetStd, RegulationType FROM Regulation WHERE  RegulatorSys LIKE '{}' OR RegulatorStd LIKE '{}'".format(regulator,regulator)
        targets = pd.read_sql_query(query,conn)
        conn.close()
        return targets
    return None
