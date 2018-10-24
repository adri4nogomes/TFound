from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from PositionMatrix import PositionMatrix as PM
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from Bio.Restriction import *
import PositionMatrix as pm
from Bio.Seq import Seq
from Bio import SeqIO
from numba import jit
import pandas as pd
import numpy as np
import random
import os

global base, tab, S2A, A2S
base = ["A","T","G","C"]
tab = str.maketrans("atgcATGC", "tacgTACG")
S2A = str.maketrans("atgcnxATGCNX", "012344012344")
A2S = str.maketrans("012344", "ATGCNX")

def getScorePostion(seq, tf, threshold, normalized=True, outOfRange=False, reversePosition=False):
    scores = []
    aux = pm.getScore(np.array(DNA2Array(getAnti(seq) if reversePosition else seq), dtype=np.int64), tf.pm, outOfRange=outOfRange, reversePosition=reversePosition)
    scores.extend(aux)
    scores = np.array(scores)
    if(normalized):
        TotalScore=[tf.getWorstMotif(first=True)[1], tf.getBestMotif(first=True)[1]]
        scores[:,0] = [-(score/TotalScore[0]) if score<0 else score/TotalScore[1] for score in scores[:,0]]
    return scores[np.where(scores[:,0]>threshold)]

#Pegar genoma (em dict)
def getGenome(file):
    with open(file, "rU") as handle:
        return dict([record.id, str(record.seq)] for record in SeqIO.parse(handle, "fasta"))

def getAnnotation(file):
    #Abrir gff
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_table(file, low_memory=False, comment='#', header=None, names=columns)
    return df

def getSequences(genome, gff, upstream, downstream=0):
    down = downstream
    gff = gff[gff.type == 'gene']
    Names = [gff.iloc[i,8].split(";")[2].split("=")[1] for i in range(len(gff))]
    columns = ["gene", "sequence", "upstream", "downstream"]
    table=[]
    for i in range(len(gff)):
        chromosome = gff.iloc[i,0]
        strand = gff.iloc[i,6]
        start = gff.iloc[i,3] if strand=="+" else gff.iloc[i,4]
        down = down if(downstream>=0) else abs(gff.iloc[i,4]-gff.iloc[i,3]+1)
        index = start-1
        inf,sup = [index-down+1,index+upstream+1] if strand=="-" else [index-upstream,index+down]
        sequence = genome[chromosome][max(0,inf):min(sup,len(genome[chromosome])-1)]
        sequence = sequence if strand=="+" else getAnti(sequence)
        u=min(sup,len(genome[chromosome])-1)-index-1 if strand=="-" else index-max(0,inf)
        d=index-max(0,inf)+1 if strand=="-" else min(sup,len(genome[chromosome])-1)-index
        table.append([Names[i], sequence, u, d])
    df = pd.DataFrame(table, columns=columns)
    df = df.sort_values(by=['gene'], ascending = True)
    return df

def RE(file, re, upstream="", core=""):
    with open(file) as f:
        sequences = f.read().rstrip().split("\n")
    if(isinstance(re, str)):
        with open(re) as f:
            aux = f.read().rstrip().split("\n")
        re=aux
        print(re)
    amb = IUPACAmbiguousDNA()
    sense = [Seq(upstream+sequence+core, amb) for sequence in sequences]
    anti = [Seq(getAnti(upstream+sequence+core), amb) for sequence in sequences]
    k=0
    for i in range(len(sequences)):
        for j in re:
            if(eval(j).search(sense[i])!=[]):
                print(sense[i] + " sub: "+ sequences[i]+": "+str(i)+" "+j+" (Sense)")
                print(eval(j).search(sense[i]))
                k +=1
            elif(eval(j).search(anti[i])!=[]):
                print(anti[i] + " sub: "+ sequences[i]+": "+str(i)+" "+j+" (Anti)")
                print(eval(j).search(anti[i]))
                k +=1
    print(str(k)+" sequências bateram com as enzimas de restrição.")
def getAnti(sequence):
    return sequence.translate(tab)[::-1]

def generateRandomSequence(length):
    return ''.join(random.choices(base, k=length))

def generateSemiRandomSequence(values):
    #completar os motifs pequenos
    dif = values[1] - len(values[0])
    if(dif>0):
        r = random.randint(0, dif)
        result = ''.join(random.choices(base, k=r)) + values[0] + ''.join(random.choices(base, k=dif-r))
        return result
    return values[0]

# gerar primeira população aleatória ou com base nos tfs
def generatePopulation(number, length, tfs=None):
    population = []
    if(tfs==None):
        population = [generateRandomSequence(length) for i in range(number)]
    else:
        L = list()
        for tf in tfs:
            #melhorar divisão
            L.extend([[tf.getBestMotif()[0],length]] * (number//len(tfs)))
        while(number-len(L)>0):
            L.append([tfs[len(tfs)-1].getBestMotif()[0],length])
        population = [generateSemiRandomSequence(l) for l in L]
    return population

def getOutputName(output):
    # Não sobrescrever arquivos existentes
    extension = output.split(".")[-1:][0]
    name = output[:-(len(extension)+1)]
    i = 0
    while os.path.exists(output):
        output = "{}[{}].{}".format(name, i, extension)
        i += 1

    return output

def get_cmap(n, name='hsv'):
    return plt.cm.get_cmap(name, n)

def PlotAllScores(sequence, tfs=[], outOfRange=False, principalOnly=False, margin=0, output=None, save=True, normalized=True, paralogs=[], maximize=[], name=None):
    db = "TF.db"
    if(os.path.isfile(db)):
        tfs_stack = [np.vstack((tf.pm for tf in tfs)),[len(tf.pm) for tf in tfs]]#Pilha de tfs e tamanhos

        sense = pm.getScores(np.array([DNA2Array(sequence)], dtype=np.int64), tfs_stack[0], tfs_stack[1], outOfRange=outOfRange, reverse = False)[0]
        anti = pm.getScores(np.array([DNA2Array(sequence)], dtype=np.int64), tfs_stack[0], tfs_stack[1], outOfRange=outOfRange, reverse = True)[0]

        if(normalized):
            TotalScore=[[tf.getWorstMotif(first=True)[1], tf.getBestMotif(first=True)[1]] for tf in tfs]
            for i in range(len(tfs)):
                sense[i] = sense[i]/TotalScore[i][1] if sense[i]>0 else -(sense[i]/TotalScore[i][0])
                anti[i] = anti[i]/TotalScore[i][1] if anti[i]>0 else -(anti[i]/TotalScore[i][0])

        # Gerar ponto no gráfico
        Min = float("inf")
        points = pd.DataFrame(columns=["TF", "Sense", "Anti-sense", "Color", "ExpertConfidence"])
        for i in reversed(range(len(tfs))):
            tam = 20.0 if(tfs[i].EC=="High") else (15.0 if(tfs[i].EC=="Medium") else (10.0 if(tfs[i].EC=="Low") else 5))

            if(i<len(maximize)): #problema
                if(maximize[i]==1):
                    Min = sense[i] if(sense[i]<Min) else Min
                elif(maximize[i]==0):
                    Min = anti[i] if(anti[i]<Min) else Min
                else:
                    Min = Min if(anti[i]>Min and sense[i]>Min) else (sense[i] if(sense[i]<anti[i]) else anti[i])
                points.loc[len(points)] = [tfs[i].StandardName, sense[i], anti[i], "R", tam]

            elif(tfs[i].StandardName != None and (tfs[i].StandardName in paralogs or tfs[i].SystematicName in paralogs or (tfs[i].MotifID!=None and tfs[i].MotifID in paralogs))):
                points.loc[len(points)] = [tfs[i].StandardName, sense[i], anti[i], "G", tam]
            else:
                points.loc[len(points)] = [tfs[i].StandardName, sense[i], anti[i], "B", tam] if(principalOnly==False) else ["", sense[i], anti[i], "B", tam]

        fig = Figure(figsize=(6,6))
        a = fig.add_subplot(111)
        a.scatter(points.loc[:,"Sense"], points.loc[:,"Anti-sense"], marker='.', c=points.loc[:,"Color"], cmap=plt.get_cmap('Spectral'), s=points.loc[:,"ExpertConfidence"])

        # Nomear pontos
        for label, px, py in zip(points.loc[:,"TF"], points.loc[:,"Sense"], points.loc[:,"Anti-sense"]):
            a.annotate(label, xy=(px, py), xytext=(0, 0), textcoords='offset points', ha='right', va='baseline', size=10, alpha=1)

        a.axhline(0, color='k', lw=0.5)
        a.axvline(0, color='k', lw=0.5)

        if(len(maximize)>0):
            a.axhline(Min, color='R', lw=0.5, ls="-")
            a.axvline(Min, color='R', lw=0.5, ls="-")

        if(margin>0):
            a.axhline(Min-margin, color='B', lw=0.5)
            a.axvline(Min-margin, color='B', lw=0.5)

        a.set_title ("Score distribution of "+name if name!=None else Array2DNA(sequence), fontsize=15, fontweight=0, color='orange')
        a.set_ylabel("Sequence' Score", fontsize=14)
        a.set_xlabel("Sequence Score", fontsize=14)

        #Não funciona o save
        if(save==True):
            output = "imgs/Scores_"+"_".join(t.upper() for t in targets)+"_"+str(margin)+".svg" if output==None else output
            output = getOutputName(output)
            os.makedirs(os.path.dirname(output), exist_ok=True)
            #Salvar plot
            a.savefig(output, bbox_inches='tight', format='svg')
            #plt.savefig(output, bbox_inches='tight', format='png', dpi=300)
            #Limpar plot
            a.gcf().clear()
            #plt.show()
        else:
            return fig
    else:
        print("ERROR: Invalid entry.")
        return None

def Array2DNA(array):
    return ''.join(str(x) for x in array).translate(A2S)

def DNA2Array(sequence):
    return list(map(int, sequence.translate(S2A)))

@jit(nopython=True)
def getFitness(sequences, pms, lpm, weight,  maximize, bestScore, worstScore, normalized=False, outOfRange=False, margin=0.0, reverse = False, maxRep=5):
    llpm = len(lpm)
    lmaximize = len(maximize)
    lsequences = len(sequences)

    fitness = np.zeros((lsequences,llpm+1), dtype=np.float64)
    S = np.zeros(lsequences, dtype=np.float64)
    D = np.zeros(lsequences, dtype=np.float64)

    Out = 0.0
    if(normalized):
        Out = 0.1
    else:
        Out = 1.0

    # Cálculo dos scores
    sense = pm.getScores(sequences, pms, lpm, outOfRange=outOfRange, reverse = False)
    antisense = pm.getScores(sequences, pms, lpm, outOfRange=outOfRange, reverse=True)

    # Normalização
    if(normalized):
        normalization(sense, bestScore, worstScore, lpm)
        normalization(antisense, bestScore, worstScore, lpm)

    for s in range(lsequences):
        if(maxRep>0 and not limitRepetitions(sequences[s], maxRep)):
            if(normalized):
                fitness[s,0] -= 2*(len(lpm)-len(maximize))
            else:
                fitness[s,0] -= 30*(len(lpm)-len(maximize))
        else:
            TFmin = 0.0
            TFmax = 0.0
            # Maior score por orientação
            if(lmaximize==0):
                if(normalized):
                    D[s] = 1.0
                    TFmin = 1.0
                else:
                    D[s] = np.amax(fitness[s,:])
                    TFmin = np.amax(fitness[s,:])
                S[s] = 1.0
            limit = -np.inf
            for l in range(llpm):
                if(l<lmaximize):
                    if(maximize[l]==1):
                        fitness[s,l+1] = sense[s,l]
                    elif(maximize[l]==0):
                        fitness[s,l+1] = antisense[s,l]
                    else:
                        fitness[s,l+1] = max(sense[s,l], antisense[s,l])
                    D[s] += fitness[s, l+1]*weight[l]
                    S[s] += weight[l]
                else:
                    if(l==lmaximize):
                        TFmin = np.amin(fitness[s, 1:lmaximize+1])
                        if(margin>0.0):
                            limit = TFmin-margin
                        else:
                            limit = TFmin

                    if(reverse==True):
                        fitness[s,l+1] = max(sense[s,l], antisense[s,l])
                    else:
                        fitness[s,l+1] = sense[s,l]

                    if(fitness[s, l+1]>limit):
                        D[s] -= np.sqrt(np.power(fitness[s, l+1]*weight[l]-limit, 2))
                        S[s] += weight[l]
                        if(fitness[s, l+1]>=TFmin):
                            fitness[s,0] -= Out
                        else:
                            fitness[s,0] -= Out/2
            TFmax = np.amax(fitness[s, lmaximize+1:])

            # Cálculo de fitness
            fitness[s,0] += D[s]/S[s]
            if(llpm>lmaximize and (TFmin-TFmax)>0.0 and margin>=0.0):
                fitness[s,0] += TFmin-TFmax

    return fitness

@jit(nopython=True)
def limitRepetitions(sequence, maxRep):
    k = np.zeros(4, dtype=np.int64)
    Max=0
    for i in range(len(sequence)):
        for j in range(4):
            if(sequence[i]==j):
                k[j] += 1
            else:
                k[j] = 0
            if(k[j]>Max):
                Max = k[j]
        if(Max>maxRep):
            return False
    return True

@jit(nopython=True)
def normalization(score, bestScore, worstScore, lpm):
    for s in range(len(score)):
        for l in range(len(lpm)):
            if(score[s,l]>0):
                score[s,l] = score[s,l]/bestScore[l]
            else:
                score[s,l] = -(score[s,l]/worstScore[l])

@jit(nopython=True)
def getMutation(sequence, mutation_rate):
    if(np.random.uniform(0.0, 1.0)<=mutation_rate):
        p = np.random.randint(len(sequence))
        b = sequence[p]
        while b==sequence[p]:
            b = np.random.randint(3)
        sequence[p] = b
    return sequence

# Retornar o cruzamento de múltiplos pais ou os próprios pais
@jit(nopython=True)
def getCrossingover(parents, crossingover_rate, individual_length):
    if(np.random.uniform(0.0, 1.0)<=crossingover_rate):
        sons = np.zeros((2,individual_length), dtype=np.int64)
        for i in range(individual_length):
            if(np.random.randint(0,1)==0):
                sons[0,i] = parents[0,i]
                sons[1,i] = parents[1,i]
            else:
                sons[0,i] = parents[1,i]
                sons[1,i] = parents[0,i]
        return sons
    else:
        return parents

@jit(nopython=True)
def getReproduction(population, size, individual_length, crossingover_rate, mutation_rate):
    npopulation = np.zeros((size, individual_length),dtype=np.int64)
    #Crossingover+Mutation
    for i in range(size-1):
        parents = np.zeros((2,individual_length), dtype=np.int64)
        for j in range(2):
            parents[j] = population[np.random.randint(len(population))]
        child = getCrossingover(parents, crossingover_rate, individual_length)
        npopulation[i,:] = getMutation(child[0], mutation_rate)
        i+=1
        npopulation[i,:] = getMutation(child[1], mutation_rate)
    return npopulation
