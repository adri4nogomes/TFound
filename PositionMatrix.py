from numba import jit, int64, int32, int8, float32, float64, boolean
from io import StringIO
import numba as nb
import pandas as pd
import numpy as np
import math
import time
import os
#cuda.select_device(0)#usar gpu extra

global base
base = ['A', 'T', 'G', 'C']

#@jit(nopython=True, parallel=True, fastmath=True)
@jit(nopython=True)
def getScore(sequence, pm, outOfRange=False, reversePosition=False):
    P=0
    newsequence = np.copy(sequence)
    #newsequence = np.concatenate((sequence, sequence, sequence))
    if(outOfRange==True or len(newsequence)<pm.shape[0]):
        P = 5-pm.shape[0]
        x = np.full(-P, 4, dtype=np.int64)
        z = np.hstack((x, newsequence, x))
        newsequence = z

    dif = len(newsequence)-pm.shape[0]+1
    #Max = np.full((dif,2), -np.inf, dtype=np.float64)
    Max = np.full((dif,2), 0, dtype=np.float64)

    for i in range(dif):
        for j in range(pm.shape[0]):
            Max[i,0] += pm[j, newsequence[j+i]]
        if(reversePosition):
            Max[i,1]=len(sequence)-(i+P)-1
        else:
            Max[i,1]=i+P

    Max[:,1] = Max[np.argsort(Max[:,0]),1]
    Max[:,0] = np.sort(Max[:,0])

    #p = np.argmax(Max[:,0]) # argumento com maior score

    return Max
@jit(nopython=True)
def _getScore(sequence, pm, outOfRange=False, reversePosition=False):
    P=0
    newsequence = np.copy(sequence)
    #newsequence = np.concatenate((sequence, sequence, sequence))
    if(outOfRange==True or len(newsequence)<pm.shape[0]):
        P = 5-pm.shape[0]
        x = np.full(-P, 4, dtype=np.int64)
        z = np.hstack((x, newsequence, x))
        newsequence = z

    dif = len(newsequence)-pm.shape[0]+1
    Max = -np.inf
    for i in range(dif):
        aux = 0
        for j in range(pm.shape[0]):
            aux += pm[j, newsequence[j+i]]
        if(aux>Max):
            Max=aux
            if(reversePosition):
                p=len(sequence)-(i+P)-1
            else:
                p=i+P

    return Max,p

@jit(nopython=True)
def getScores(sequences, pms, lpm, outOfRange=False, reverse=False):
    scores = np.full((len(sequences),len(lpm)), -np.inf, dtype=np.float64)
    for i in range(len(sequences)):
        ini = 0
        end = 0
        if(reverse==True):
            antisense = np.copy(sequences[i])[::-1]
            for l in range(len(antisense)):
                if(antisense[l]==0):
                    antisense[l] = 1
                elif(antisense[l]==1):
                    antisense[l] = 0
                elif(antisense[l]==2):
                    antisense[l] = 3
                elif(antisense[l]==3):
                    antisense[l] = 2
            for j in range(len(lpm)):
                end += lpm[j]
                scores[i,j] = _getScore(antisense, pms[ini:end,:], outOfRange)[0]
                ini = end
        else:
            for j in range(len(lpm)):
                end += lpm[j]
                scores[i,j] = _getScore(sequences[i], pms[ini:end,:], outOfRange)[0]
                ini = end

    return scores

def PFMtoPPM(pm):
    row_sums = pm.sum(axis=1)
    pm = pm / row_sums[:,np.newaxis]
    return pm

#Transforma arquivos de sequências, strings de sequências ou listas de sequências em PFM
def seqsToPFM(sequences, sep="\n"):
    #Transformar string em lista
    try:
        if isinstance(sequences, str):
            #Abrir arquivo
            if(sequences.count('.')>=1):
                with open(sequences) as f:
                    sequences = f.read().rstrip().split(sep)
            #Quebrar string
            else:
                sequences = sequences.split(sep)
        sequences = [sequence.strip() for sequence in sequences]
        l = len(sequences[0])

        if isinstance(sequences, list) and all(isinstance(sequence, str) and len(sequence)==l for sequence in sequences):
            pfm = np.zeros(shape=(l,4))
            for sequence in sequences:
                for i in range(len(sequence)):
                    pfm[i, base.index(sequence[i])] += 1
            return pfm
        else:
            print("ERROR: Invalid entry.")
    except Exception as e:
        raise e

class PositionMatrix (object):
    """Armazena um fator de transcrição e todas suas propriedades"""
    def __init__(self, pm, StandardName=None, SystematicName=None, MotifID=None, SubMotif=None, sep="\t", type="pwm", EC=None, TotalScore=None):
        try:
            self.StandardName = StandardName
            self.SystematicName = SystematicName
            self.MotifID = MotifID
            self.SubMotif = SubMotif
            self.EC = EC
            self.type = type
            self.TotalScore=TotalScore

            if isinstance(pm, np.ndarray):
                self.pm = pm
            elif isinstance(pm, pd.DataFrame):
                self.pm = pm.as_matrix()
            elif isinstance(pm, str):
                self.pm = self.String2PM(pm, sep)
            else:
                print("ERRO: Parâmetro inválido.")
        except Exception as e:
            raise e

    def getMaxScore(self, sequence, outOfRange=False):
        #retornar o maior score para uma sequência e a posição inicial do encaixe
        return list(getScore(np.array(sequence, dtype=np.int64), self.pm, outOfRange))

    def PFMtoPPM(self):
        self.pm = PFMtoPPM(self.pm)
        self.type = "ppm"

    def getBestMotif(self, first=False):
        l = [[],0]
        self.getMotifR(l, best=True, first=first)
        return l

    def getMotifR(self, *l, best=None, i=0, L=[], S=0, first=False):
        if(i==self.pm.__len__()):
            if(len(l[0][0])==0):
                l[0][1] = S
            l[0][0].append(''.join(L))
        else:
            target =  np.amax(self.pm[i]) if best else np.amin(self.pm[i])
            for b in range(len(base)):
                found=False
                if(self.pm[i, b]==target):
                    x = L.copy()
                    x.append(base[b])
                    self.getMotifR(*l, best=best, i=(i+1), L=x, S=(S+target), first=first)
                    found=True
                if(found and first): break

    def getWorstMotif(self, first=False):
        l = [[],0]
        self.getMotifR(l, best=False, first=first)
        return l

    def __len__(self):
        #retornar tamanho do pwm
        return self.pm.shape[0]

    def __str__(self):
        df = pd.DataFrame(data=self.pm[:,:-1],index=range(self.__len__()),columns=base)
        return "\n"+str(self.name())+"\n" + str(df)

    def __repr__(self):
        return self.__str__()

    def name(self):
        if(self.MotifID!=None):
            return self.MotifID
        elif(self.StandardName!=None):
            return self.StandardName
        elif(self.SystematicName!=None):
            return self.SystematicName
        else:
            return None

    def String2PM(self, string, sep):
        # Conversão da string no formato do site para dataframe
        try:
            if isinstance(string, str):
                if(string.count('\n')<=1):
                    self.type = string.split(".")[-1]
                    matrix = np.genfromtxt(string, delimiter=sep)
                else:
                    matrix = np.genfromtxt(StringIO(string), delimiter=sep)
                matrix = matrix.transpose()
                matrix = np.delete(matrix, (0), axis=0)
                # Adicionar coluna do mínimo
                matrix = np.hstack((matrix, np.amin(matrix, axis=1).reshape(matrix.shape[0], 1)))
                return matrix
            else:
                print("ERROR: Invalid entry.")
        except Exception as e:
            raise e

    def SaveFile(self, output=None):
        try:
            # Gerar nome
            if(output == None):
                output = "pm/" + self.name + "." + self.type

            # Não sobrescrever arquivos existentes
            extension = output.split(".")[-1:][0]
            name = output[:-(len(extension)+1)]
            i = 0
            while os.path.exists(output):
                output = "{}[{}].{}".format(name, i, extension)
                i += 1
            # Criar diretórios
            os.makedirs(os.path.dirname(output), exist_ok=True)

            with open(output, 'w') as file:
                for b in base:
                    line = self.pm[:,base.index(b)].tolist()
                    file.write(b),
                    for l in line:
                        file.write("\t" + str(l)),
                    file.write("\n")
        except Exception as e:
            raise e
