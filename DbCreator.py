from urllib.request import urlopen
from urllib.parse import urlencode
from bs4 import BeautifulSoup
import multiprocessing as mp
from io import StringIO
import pandas as pd
import sqlite3
import tarfile
import shutil
import time
import json
import re
import os

#Download de arquivo
def Download(url, folder="."):
    file_name = os.path.join(folder, url.split('/')[-1]) #Gerar nome do arquivo
    if(os.path.isfile(file_name)): os.remove(file_name) #Excluir arquivo existente
    u = urlopen(url) #Carregar página
    f = open(file_name, 'wb') #Criar arquivo de bytes
    #Download
    while True:
        buffer = u.read(8192)
        if not buffer: break
        f.write(buffer)
    f.close()
    print("Download completed.")
    return file_name

#Extrai e remove arquivo
def ExtractTar(file_name, folder="."):
    #Descompactação
    tar = tarfile.open(file_name, "r:gz")
    tar.extractall(path=folder)
    tar.close()
    time.sleep(2)
    os.remove(file_name)
    print("Extraction completed.")

#Cria o banco de dados de arquivos dos fatores de transcrição
def DownloadDB():
    folder = "PWMS"
    if(os.path.isdir(folder)): shutil.rmtree(folder) #Excluir pasta temporária se ela já existe
    os.makedirs(folder) #Criar pasta temporária para salvar

    #Fazer download do arquivo para a pasta
    #ExtractTar(Download("http://yetfasco.ccbr.utoronto.ca/1.02/Downloads/Expert_PWMs.tar.gz", folder=folder), folder) #Curated
    ExtractTar(Download("http://yetfasco.ccbr.utoronto.ca/1.02/Downloads/All_PWMs.tar.gz", folder=folder), folder) #All

    return SearchFiles(folder, "YeTFaSCo")

def SearchFiles(folder, db):
    #Vasculhar pasta em busca de arquivos
    Exts = []
    for (dirpath, dirnames, filenames) in os.walk(folder):
        if(filenames!=[]):
            for file in filenames:
                aux = file.split(".")
                name, ext = ['.'.join(aux[:-1]),aux[-1]] if len(aux)>1 else [file,None] #Pegar nome arquivo e extensão
                Exts.append([os.path.join(os.path.abspath(dirpath), file), ext, name, db])
    return Exts

def files2DB(db, Exts, New=False):
    #Criação do BD
    if(os.path.isfile(db) and New): os.remove(db)
    conn = sqlite3.connect(db)
    c = conn.cursor()

    if(New):
        #Criação das tabelas
        c.execute('''CREATE TABLE Files(absPath TEXT NOT NULL UNIQUE, Extension TEXT,
            Name TEXT NOT NULL, Database TEXT NOT NULL, PRIMARY KEY (absPath))''') #Criar tabela de extensões

    #Inserir dados
    c.executemany('INSERT INTO Files(absPath, Extension, Name, Database) VALUES(?,?,?,?)', Exts)

    conn.commit() #Salvar alterações
    conn.close() #Fechar conexão
    print("Database created.")
    return [Ext[2] for Ext in Exts]

#Completa informações de fatores de transcrição de leveduras com auxílio do SGD e do YeTFaSCo
def YeastInfo(db, names):
    if(os.path.isfile(db)):
        tfs = set()
        for name in names:
            try:
                if(len(name.split("_")[0].split("-")[1])==1): raise Exception()
            except: tfs.add(name.split("_")[0])

        #Buscar info no YeTFaSCo
        tfFiles = YeTFaSCo(names)

        p = mp.Pool(mp.cpu_count()*4)
        #Buscar info no SGD
        tfSummary = p.map_async(summarySGD, tfs).get()
        tfRegulation = p.map_async(regulationSGD, tfs).get()
        p.close()
        p.join()

        tfSummary = [tf for tf in tfSummary if tf!=None ]
        TFs = [tuple(tf[:-1]) for tf in tfSummary]
        Paralogs = set(frozenset([tf[0], paralog]) for tf in tfSummary for paralog in tf[-1])
        Paralogs = [tuple(item) for item in Paralogs]

        tfRegulation = [tf for tf in tfRegulation if tf!=None ]
        tfRegulation = [item for sublist in tfRegulation for item in sublist]

        #Criar tabelas
        conn = sqlite3.connect(db)
        c = conn.cursor()

        #Criação das tabelas
        c.execute('''CREATE TABLE TFs(SystematicName TEXT NOT NULL, StandardName TEXT,
            SGDId TEXT UNIQUE, Sense INTEGER NOT NULL, Antisense INTEGER NOT NULL,
            PRIMARY KEY (SystematicName,StandardName))''')
        c.execute('''CREATE TABLE FileTF(SystematicName TEXT NOT NULL, fileName TEXT NOT NULL,
            MotifID INTEGER, SubMotif INTEGER, TotalScore REAL,
            DBS TEXT, ExpertConfidence TEXT, MethodDescription TEXT, PMID INTEGER,
            Dubious BOOLEAN, Dimer BOOLEAN NOT NULL, FOREIGN KEY(fileName) REFERENCES Files(Name),
            PRIMARY KEY (SystematicName, fileName))''')
        c.execute('''CREATE TABLE Paralogs(ParalogA TEXT NOT NULL, ParalogB TEXT NOT NULL,
            FOREIGN KEY(ParalogA) REFERENCES TFs(SystematicName),
            PRIMARY KEY (ParalogA, ParalogB))''')
        c.execute('''CREATE TABLE Regulation(RegulatorSys TEXT NOT NULL, RegulatorStd TEXT NOT NULL,
            TargetSys TEXT NOT NULL, TargetStd TEXT NOT NULL, RegulationType TEXT NOT NULL, Evidence TEXT, HappensDuring TEXT,
            Strain TEXT NOT NULL, FOREIGN KEY(TargetSys) REFERENCES TFs(SystematicName),
            PRIMARY KEY (RegulatorSys, TargetSys))''')

        #Inserir dados
        c.executemany('INSERT INTO  TFs(SystematicName, StandardName, SGDId, Sense, Antisense) VALUES(?,?,?,?,?)', TFs)
        c.executemany('INSERT INTO Paralogs(ParalogA, ParalogB) VALUES(?,?)', Paralogs)
        c.executemany('''INSERT INTO  FileTF(SystematicName, fileName, MotifID, SubMotif, TotalScore, DBS, ExpertConfidence,
            MethodDescription, PMID, Dubious, Dimer) VALUES(?,?,?,?,?,?,?,?,?,?,?)''', tfFiles)
        c.executemany('''INSERT OR IGNORE INTO Regulation(RegulatorSys, RegulatorStd, TargetSys, TargetStd, RegulationType, Evidence, HappensDuring, Strain)
            VALUES(?,?,?,?,?,?,?,?)''', tfRegulation)

        conn.commit() #Salvar alterações
        conn.close() #Fechar conexão
    else:
        print("ERROR: DB not found.")

def EmptyInfo(db, names):
    if(os.path.isfile(db)):
        #Gerar dados para inserção
        TFs = []
        tfFiles = []

        MotifID=None
        SubMotif=None
        TotalScore=None
        DBS=None
        ExpertConfidence=None
        MethodDescription=None
        PMID=None
        Dubious=None
        Dimer=False
        SGDId=None
        sense=0
        antisense=0

        for name in names:
            SystematicName=name
            StandardName=name
            fileName=name
            TFs.append((SystematicName, StandardName, SGDId, sense, antisense))
            tfFiles.append((SystematicName, fileName,MotifID,SubMotif,TotalScore,DBS,ExpertConfidence,MethodDescription,PMID,Dubious,Dimer))

        #Criar tabelas
        conn = sqlite3.connect(db)
        c = conn.cursor()

        #Inserir dados
        c.executemany('INSERT INTO  TFs(SystematicName, StandardName, SGDId, Sense, Antisense) VALUES(?,?,?,?,?)', TFs)
        c.executemany('''INSERT INTO  FileTF(SystematicName, fileName, MotifID, SubMotif, TotalScore, DBS, ExpertConfidence,
            MethodDescription, PMID, Dubious, Dimer) VALUES(?,?,?,?,?,?,?,?,?,?,?)''', tfFiles)

        conn.commit() #Salvar alterações
        conn.close() #Fechar conexão
    else:
        print("ERROR: DB not found.")

def YeTFaSCo(names):
    try:
        html = urlopen("http://yetfasco.ccbr.utoronto.ca/MotViewLong.php?PME_sys_export=csv")
    except:
        return None
    else:
        #Converter para pandas
        csv = []
        for line in html.read().decode('utf-8').split("\n"): csv.append(line[1:].replace('"', ''))
        table = pd.read_csv(StringIO('\n'.join(csv)), sep=';')

        result = []
        for name in names:
            #Pegar nome e id
            SystematicName, MotifID = name.split("_")
            try: MotifID, SubMotif = MotifID.split(".")
            except: SubMotif = 0

            item=table.loc[(table['Motif ID']==int(MotifID)) & (table['Sub Motif']==int(SubMotif))]
            #Buscar dados por motif id e submotif
            StandardName=item['Systematic Name'].values[0]
            TotalScore=item['Total Score'].values[0]
            DBS=item['DBDs'].values[0]
            ExpertConfidence=item['Expert Confidence'].values[0]
            MethodDescription=item['Method Description'].values[0]
            PMID=int(item['Reference'].values[0])
            Dubious=True if item['Dubious?'].values[0]=="Dubious" else False

            try:
                aux = SystematicName.split("-")
                if(len(aux[1])==1): raise Exception()
                SystematicName = aux if aux[-1]!="dimer" else aux[:-1]
                result.extend([(SN.upper(),name,MotifID,SubMotif,TotalScore,DBS,ExpertConfidence,MethodDescription,PMID,Dubious,True) for SN in SystematicName])
            except:
                result.append((SystematicName,name,MotifID,SubMotif,TotalScore,DBS,ExpertConfidence,MethodDescription,PMID,Dubious,False))
        return result

def summarySGD(name):
    while True:
        try:
            html = urlopen("https://www.yeastgenome.org/webservice/locus/{}/".format(name))
        except:
            print("{} not found in the SGD.".format(name))
            return None
        else:
            res = json.loads(html.read().decode('utf-8')) # carregando da string

            StandardName = res['display_name']
            SystematicName = res['format_name']
            SGDId = res['sgdid'] #Saccharomyces Genome Database ID
            Paralogs = [paralog['child']['format_name'] for paralog in res['paralogs']]
            BP = [bp['term']['display_name'] for bp in res['go_overview']['manual_biological_process_terms']]

            sense = 0
            antisense = 0
            for bp in BP:
                if(re.search("positive regulation of antisense(.*)RNA", bp)):
                    antisense += 1 if antisense==0 or antisense==2 else 0
                elif(re.search("negative regulation of antisense(.*)RNA", bp)):
                    antisense += 2 if antisense==0 or antisense==1 else 0
                elif(re.search("positive regulation of transcription", bp) or re.search("^transcription", bp)
                    or re.search("activation of transcription", bp) or re.search("RNA(.*)preinitiation complex assembly", bp)):
                    sense += 1 if sense==0 or sense==2 else 0
                elif(re.search("negative regulation(.*)RNA", bp)):
                    sense += 2 if sense==0 or sense==1 else 0

            return [SystematicName, StandardName, SGDId, sense, antisense, Paralogs]

def regulationSGD(name):
    while True:
        try:
            html = urlopen("https://www.yeastgenome.org/webservice/locus/{}/regulation_details".format(name))
        except:
            print("{} not found in the SGD.".format(name))
            return None
        else:
            res = json.loads(html.read().decode('utf-8')) # carregando da string

            Regulation = []
            for reg in res:
                RegulationType = reg['regulation_type']
                Evidence = reg['evidence']['display_name']
                RegulatorSys = reg['locus1']['format_name']
                RegulatorStd = reg['locus1']['display_name']
                TargetSys = reg['locus2']['format_name']
                TargetStd = reg['locus2']['display_name']
                HappensDuring = reg['happens_during']
                Strain = reg['strain']['display_name']
                Regulation.append((RegulatorSys, RegulatorStd, TargetSys, TargetStd, RegulationType, Evidence, HappensDuring, Strain))

            return Regulation

def main():
    db = 'TF.db'
    names = files2DB(db, DownloadDB(), New=True)
    YeastInfo(db, names)

if __name__ == '__main__':
    main()

def YeTFaSCoOLD(name):
    #http://yetfasco.ccbr.utoronto.ca/MotViewLong.php?PME_sys_export=csv

    #Pegar nome e id
    SystematicName, MotifID = name.split("_")
    try:
        MotifID, SubMotif = MotifID.split(".")
    except:
        SubMotif = 0

    #post_args = urlencode({'PME_sys_qf1':SystematicName,'PME_sys_qf2':MotifID, 'PME_sys_qf3':SubMotif}).encode("utf-8")
    post_args = urlencode({'PME_sys_qf2':MotifID, 'PME_sys_qf3':SubMotif}).encode("utf-8")

    while True:
        try:
            html = urlopen("http://yetfasco.ccbr.utoronto.ca/MotViewLong.php", post_args)
        except:
            return None
        else:
            res = BeautifulSoup(html.read(),"html5lib");
            aux = res.findAll("td", {"class": "pme-cell-0"})

            StandardName = aux[0].getText() if(aux[0].getText()!="\xa0") else None
            TotalScore = aux[5].getText()
            DBS = aux[6].getText() if(aux[6].getText()!="\xa0") else None
            ExpertConfidence = aux[7].a.getText() if aux[7].a!=None else None
            MethodDescription = aux[8].getText()
            PMID = aux[9].getText()
            Dubious = True if aux[10].getText()=="Dubious" else False
            try:
                aux = SystematicName.split("-")
                if(len(aux[1])==1):
                    raise Exception()
                SystematicName = aux if aux[-1]!="dimer" else aux[:-1]
                return [(SN,name,MotifID,SubMotif,TotalScore,DBS,ExpertConfidence,MethodDescription,PMID,Dubious,True) for SN in SystematicName]
            except:
                return [(SystematicName,name,MotifID,SubMotif,TotalScore,DBS,ExpertConfidence,MethodDescription,PMID,Dubious,False)]
