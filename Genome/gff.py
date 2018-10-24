import pandas as pd
#import GeneticAlgorithm as ga
from Bio import SeqIO
#import juju as j
import numpy as np

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
        #dar ga.getAnti na seq "-"
        u=min(sup,len(genome[chromosome])-1)-index-1 if strand=="-" else index-max(0,inf)
        d=index-max(0,inf)+1 if strand=="-" else min(sup,len(genome[chromosome])-1)-index
        table.append([Names[i], sequence, u, d])
    df = pd.DataFrame(table, columns=columns)
    print(df.head())
    return df

genome_file = "GCF_000146045.2_R64_genomic.fna"
gff_file = "GCF_000146045.2_R64_genomic.gff"
#genome_file = "GCF_000007565.2_ASM756v2_genomic.fna"
#gff_file = "GCF_000007565.2_ASM756v2_genomic.gff"

genome = getGenome(genome_file)

gff = getAnnotation(gff_file)

getSequences(genome, gff, 0, -1)
'''x=1
print(gff.iloc[x])
print(spliceGenome(genome, gff.iloc[x,0], gff.iloc[x,3] if(gff.iloc[x,6]=="+") else gff.iloc[x,4], gff.iloc[x,6], 1000,3))
x=2
print(gff.iloc[x])
print(spliceGenome(genome, gff.iloc[x,0], gff.iloc[x,3] if(gff.iloc[x,6]=="+") else gff.iloc[x,4], gff.iloc[x,6], 1000,3))'''
#gi = getGeneInfo(gff)
'''print("Genes: "+str(len(gi)))

sequences = []
for i in range(len(gi)):
    item = gi.iloc[i]
    sequence = upstream(genome, item['seqid'], (item['start'] if item['strand']=="+" else item['end']), item['strand'], 1000)
    #print(">"+item['Name']+"\n"+sequence)
    sequences.append(sequence)

#pattern = "GTNCGNTNANCGNAC"
#pattern = "GTKCGMTDAWCGSAC"
#pattern = "KYKCGMWDAWCGSMC"
pattern = "GTGCGATAAACGCAC"
print(pattern)
print(ga.getAnti(pattern))
maxError = 5
motifs = j.localize(sequences, pattern, maxError=maxError, reverse=True)

for e in range(maxError+1):
    print("\nMax errors: "+str(e))

    print("Motifs: ")
    m = [motif for motif in motifs if motif[3]==e]
    print(m)
    print("Total motifs: "+str(len(m)))

    print("Genes: ")
    g = set([motif[4] for motif in m])
    print([gi['Name'].iloc[i] for i in g])
    print(g)
    print("Total genes: "+str(len(g)))'''

#Buscar por nome
#Histogram
