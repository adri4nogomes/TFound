import os
import csv
import matplotlib
import DbUtils as dbu
import pandas as pd
try:
    import Tkinter as tk
except:
    import tkinter as tk
try:
    import ttk
except:
    from tkinter import ttk
from tkinter import *
from Bio import SeqIO
import DbCreator as dbc
matplotlib.use("TkAgg")
import matplotlib as mpl
import canvasSequence as cs
import canvasPosition as cp
import Utils
from tkinter import messagebox
from tkinter import filedialog
from matplotlib.figure import Figure
from tkinter.scrolledtext import ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import platform
#from win32api import GetMonitorInfo, MonitorFromPoint

#entrada: bd pronto ou gerar bd com nome e adicionar
class InterfaceView():
    def __init__(self):
        self.OS = platform.system()
        self.master = Tk()
        self.master.title("TFound")

        #self.master.wm_iconbitmap(“endereço.ico”)
        self.master.minsize(1024, 600)

        s = ttk.Style()
        if(self.OS=="Windows"):
            s.theme_use('xpnative')
        elif(self.OS=="Linux"):
            s.theme_use('clam')
        elif(self.OS=="Darwin"):
            s.theme_use('aqua')
        # fonts for all widgets
        self.master.option_add("*Font", "arial 9")

        #Menu de opções
        self.fmenu = ttk.Frame(self.master, width=400)

        self.fcontainer = ttk.Notebook(self.fmenu)

        #Saccharomyces cerevisiae
        self.page1 = ttk.Frame(self.fcontainer)
        self.fsac = ttk.Frame(self.page1)
        self.fsacOpt = ttk.Frame(self.fsac)
        ttk.Label(self.fsacOpt, text = "Search:").pack(side=LEFT)
        ssv = StringVar()
        ssv.trace("w", lambda name, index, mode, ssv=ssv: self.searchGenes())
        self.sacSearch=ttk.Entry(self.fsacOpt, textvariable=ssv)
        self.sacSearch.pack(side=LEFT)
        ttk.Label(self.fsacOpt, text = "Regulator:").pack(side=LEFT,padx=(5,0))
        self.regulator = ttk.Combobox(self.fsacOpt, state="readonly", width=14)
        self.regulator.pack(side=LEFT)

        self.fsacOpt.pack(side=TOP, fill=X, expand=YES)
        ttk.Style().configure("Treeview.Heading", font=(None, 11))
        self.fstree = ttk.Frame(self.fsac)
        self.sactree = ttk.Treeview(self.fstree,height=7)
        self.sactree['show'] = 'headings' # remover primeira coluna vazia
        self.sactree["columns"]=("check","gene","sequence")
        self.sactree.column("check", width=1, anchor='center')
        self.sactree.column("gene", width=20, anchor='center',stretch=True)
        self.sactree.heading("check", text="☐")#" ☐ ☑ ☒ "
        self.sactree.heading("gene", text="Gene")
        self.sactree.heading("sequence", text="Sequence")
        self.sbar1= ttk.Scrollbar(self.fstree, orient='vertical', command=self.sactree.yview)
        self.sactree.configure(yscroll=self.sbar1.set)
        self.sbar1.pack(side=RIGHT, fill=Y)
        self.sactree.pack(side=TOP, fill=BOTH, expand=YES)
        self.fstree.pack(side=TOP,fil=X,expand=YES, padx=(0,12))

        self.fsave = ttk.Frame(self.fsac)
        self.lgenes = ttk.Label(self.fsave, text = "")
        self.lgenes.pack(side=LEFT)
        ttk.Button(self.fsave, text="Export", command=self.saveFile).pack(side=RIGHT, padx=2)
        self.fsave.pack(side=TOP, fil=X)

        self.fsOptions = ttk.Frame(self.fsac)
        ttk.Label(self.fsOptions, text = "Upstream (bp):").pack(side=LEFT)
        self.vint = (self.master.register(self.validate), '%P', '%s', '%S', '%W', '%V', '0', 'Int')
        self.supstream=ttk.Entry(self.fsOptions, width=6, validate="all", validatecommand=self.vint)
        self.supstream.pack(side=LEFT)
        self.supstream.insert(END, '1000')
        ttk.Label(self.fsOptions, text = "Downstream (bp):").pack(side=LEFT,padx=(5,0))
        self.sdownstream=ttk.Entry(self.fsOptions, width=6, validate="all", validatecommand=self.vint)
        self.sdownstream.pack(side=LEFT,padx=(0,5))
        self.sdownstream.insert(END, '0')
        self.sfgState = BooleanVar()
        self.sfg = ttk.Checkbutton(self.fsOptions, text = "Full gene", variable=self.sfgState, command=self.fullGene)
        self.sfg.pack(side=LEFT)
        self.fsOptions.pack(side=TOP, fill=X, expand=YES, padx=2)
        self.fsac.pack(side=TOP, fill=X)

        # Sequências
        self.page2 = ttk.Frame(self.fcontainer)
        self.fsequences = ttk.Frame(self.page2)
        self.sequence = ScrolledText(self.fsequences, height=16)
        self.sequence.insert(END, '')
        self.sequence.pack(side=TOP, fill=X, expand=YES)
        self.fsequences.pack(side=TOP, fill=X)

        #Grid de FASTA
        self.page4 = ttk.Frame(self.fcontainer)
        self.ffasta = ttk.Frame(self.page4)
        self.ffile = ttk.Frame(self.ffasta)
        ttk.Button(self.ffile, text="Load File", command=self.loadFile).pack(side=LEFT, padx=2)
        self.filePath = ttk.Label(self.ffile, text = "")
        self.filePath.pack(side=LEFT)
        self.ffile.pack(side=TOP, fill=X)
        self.fseq = ttk.Frame(self.ffasta)
        ttk.Label(self.fseq, text = "Search:").pack(side=LEFT)
        self.seqSearch=ttk.Entry(self.fseq,textvariable=ssv)
        self.seqSearch.pack(side=LEFT)
        self.lgenesF = ttk.Label(self.fseq, text = "")
        self.lgenesF.pack(side=RIGHT)

        self.fseq.pack(side=TOP, fill=X, expand=YES)
        self.ffastatree = ttk.Frame(self.ffasta)
        self.seqtree = ttk.Treeview(self.ffastatree,height=6)
        self.seqtree['show'] = 'headings' # remover primeira coluna vazia
        self.seqtree["columns"]=("check","name","sequence")
        self.seqtree.column("check", width=1, anchor='center')
        self.seqtree.column("name", width=20, anchor='center')
        self.seqtree.heading("check", text="☐")#" ☐ ☑ ☒ "
        self.seqtree.heading("name", text="Name")
        self.seqtree.heading("sequence", text="Sequence")
        self.bar1= ttk.Scrollbar(self.ffastatree, orient='vertical', command=self.seqtree.yview)
        self.seqtree.configure(yscroll=self.bar1.set)
        self.bar1.pack(side=RIGHT, fill=Y)
        self.seqtree.pack(side=TOP, fill=BOTH, expand=YES)
        self.ffastatree.pack(side=TOP,fil=X,expand=YES,padx=(0,12))

        self.fsaveF = ttk.Frame(self.ffasta)
        self.bload = ttk.Button(self.fsaveF, text="Load GFF", command=self.loadGFF)
        self.filePathGFF = ttk.Label(self.fsaveF, text = "") #Manter opções desabilitadas enquanto estiver vazio
        self.filePathGFF.pack(side=LEFT)
        ttk.Button(self.fsaveF, text="Export", command=self.saveFile).pack(side=RIGHT, padx=2)
        self.fsaveF.pack(side=TOP, fill=X)


        self.fgOptions = ttk.Frame(self.ffasta)
        ttk.Label(self.fgOptions, text = "Upstream (bp):").pack(side=LEFT)
        self.upstream=ttk.Entry(self.fgOptions, width=6, validate="all", validatecommand=self.vint)
        self.upstream.pack(side=LEFT)
        self.upstream.insert(END, '1000')
        ttk.Label(self.fgOptions, text = "Downstream (bp):").pack(side=LEFT,padx=(5,0))
        self.downstream=ttk.Entry(self.fgOptions, width=6, validate="all", validatecommand=self.vint)
        self.downstream.pack(side=LEFT,padx=(0,5))
        self.downstream.insert(END, '0')
        self.fgState = BooleanVar()
        self.fg = ttk.Checkbutton(self.fgOptions, text = "Full gene", variable=self.fgState, command=self.fullGene)
        self.fg.pack(side=LEFT)

        self.ffasta.pack(side=TOP, fill=X)

        self.fcontainer.add(self.page1, text='Saccharomyces cerevisiae')
        self.fcontainer.add(self.page2, text='Sequence')
        #self.fcontainer.add(self.page3, text='SFS')
        self.fcontainer.add(self.page4, text='FASTA/GFF')

        self.fcontainer.pack(side=TOP, fill=X)

        ttk.Separator(self.fmenu, orient="horizontal").pack(side=TOP, fill=X, pady=3)

        #DB
        self.fdb = ttk.Frame(self.fmenu)
        #Combobox DBs
        ttk.Label(self.fdb, text = "Database:").pack(side=LEFT)
        self.dbbox = ttk.Combobox(self.fdb,state="readonly",width=10)
        self.dbbox.pack(side=LEFT)
        #Carregamento de bds
        self.DBs = dbu.getDBs()
        self.dbbox['values'] = self.DBs['Database'].unique().tolist()
        self.dbbox.current(0)
        self.fdb.pack(side=TOP, fill=X)

        self.fnewdb = ttk.Frame(self.fmenu)
        ttk.Label(self.fnewdb, text = "New Database Name:").pack(side=LEFT)
        self.newdb=ttk.Entry(self.fnewdb, width=15)
        self.newdb.pack(side=LEFT)
        #ttk.Button Create New
        ttk.Button(self.fnewdb, text="Create", command=self.createDB).pack(side=LEFT)
        self.fnewdb.pack(side=TOP, fill=X)

        ttk.Separator(self.fmenu, orient="horizontal").pack(side=TOP, fill=X, pady=3)
        self.ftfs = ttk.Frame(self.fmenu)


        self.fsubthr = ttk.Frame(self.ftfs)
        ttk.Label(self.fsubthr, text = "Threshold:").pack(side=LEFT)
        self.vfloat = (self.master.register(self.validate), '%P', '%s', '%S', '%W', '%V', '0.8', 'Float')
        self.threshold = ttk.Entry(self.fsubthr, width=6, validate="all", validatecommand=self.vfloat)
        self.threshold.insert(END, '0.8')
        self.threshold.pack(side=LEFT)
        self.fsubthr.grid(row=0, column=0, sticky=E+W)

        self.normalizedState = BooleanVar()
        self.normalized = ttk.Checkbutton(self.ftfs, text = "Normalized", variable=self.normalizedState)
        self.normalized.grid(row=0, column=1, sticky=E+W)
        self.normalizedState.set(True)

        self.fsubtfs = ttk.Frame(self.ftfs)
        ttk.Label(self.fsubtfs, text = "Confidence:").pack(side=LEFT)
        self.confidence = ttk.Combobox(self.fsubtfs, state="readonly", width=8)
        self.confidence.pack(side=LEFT)
        self.confidence['values'] = self.DBs['ExpertConfidence'].loc[self.DBs['Database']==self.dbbox.get()].tolist()
        self.fsubtfs.grid(row=1, column=0, sticky=E+W)

        self.foorState = BooleanVar()
        self.outOfRange = ttk.Checkbutton(self.ftfs, text = "Out of range", variable=self.foorState)
        self.outOfRange.grid(row=1, column=1, sticky=E+W)
        self.foorState.set(True)

        ttk.Label(self.ftfs, text = "Transcription Factors:").grid(row=2, sticky=E+W)
        self.f0 = ttk.Frame(self.ftfs)
        stfv = StringVar()
        stfv.trace("w", lambda name, index, mode, stfv=stfv: self.updateSBAll())
        self.search=ttk.Entry(self.f0, textvariable=stfv, width=25)
        self.search.pack(padx=(0,12))
        self.scrollbarAll = ttk.Scrollbar(self.f0)
        self.scrollbarAll.pack(side=RIGHT, fill=Y)
        self.vAll = Variable()
        self.lAll = Listbox(self.f0, height=5, width=25, yscrollcommand=self.scrollbarAll.set, listvariable=self.vAll)
        self.lAll.pack()
        self.scrollbarAll.config(command=self.lAll.yview)
        self.f0.grid(row=3, column=0, rowspan=5, sticky=E+W)

        ttk.Label(self.ftfs, text = "Targets:").grid(row=2, column=1, sticky=E+W)
        self.f1 = ttk.Frame(self.ftfs)
        self.scrollbarTargets = ttk.Scrollbar(self.f1)
        self.scrollbarTargets.pack(side=RIGHT, fill=Y)
        self.vtargets = Variable()
        self.ltargets = Listbox(self.f1, height=6, width=25, yscrollcommand=self.scrollbarTargets.set, listvariable=self.vtargets)
        self.ltargets.pack()
        self.scrollbarTargets.config(command=self.ltargets.yview)
        self.f1.grid(row=3, column=1, rowspan=5, sticky=E+W)

        self.ftfs.pack(side=TOP)

        ttk.Button(self.ftfs, text = "Add (+) target", command= lambda: self.addTF("+")).grid(row=9, column=0, sticky=E+W)
        ttk.Button(self.ftfs, text = "Add (-) target", command= lambda: self.addTF("-")).grid(row=10,  column=0, sticky=E+W)
        ttk.Button(self.ftfs, text = "Remove target", command=self.removeTarget).grid(row=9,  column=1, sticky=E+W)
        ttk.Button(self.ftfs, text = "Reset", command=self.reset).grid(row=10,  column=1, sticky=E+W)
        ttk.Button(self.ftfs, text = "Score", command=self.Score).grid(row=11,  column=0, columnspan=2, sticky=E+W+S+N)


        self.fmenu.pack(side=LEFT, fill=Y, padx=10, pady=3)
        self.fmenu.pack_propagate(0)

        #Opções extra
        self.fcontOptions = ttk.Notebook(self.master)

        self.fsingleOpt = ttk.Frame(self.fcontOptions)
        self.fsoopt = ttk.Frame(self.fsingleOpt)
        self.allLablesState = BooleanVar()
        self.allLables = ttk.Checkbutton(self.fsoopt, text = "All ttk.Labels", variable=self.allLablesState)
        self.allLables.pack(side=LEFT)
        self.paralogState = BooleanVar()
        self.cbparalogs = ttk.Checkbutton(self.fsoopt, text = "Mark paralogs", variable=self.paralogState)
        self.cbparalogs.pack(side=LEFT)
        self.paralogState.set(True)
        self.fsoopt.pack(side=TOP, fill=X)

        self.fparalogs=ttk.Frame(self.fsingleOpt)
        self.lparalogs = Message(self.fparalogs, text = "Paralogs: {}", width=700)
        self.lparalogs.pack(side=LEFT, fill=X)
        self.fparalogs.pack(side=TOP, fill=X)

        self.fviewSingle = ttk.Frame(self.fsingleOpt)
        self.fcSequence = ttk.Frame(self.fviewSingle)
        self.fcSequence.pack(side=TOP, fill=X)
        self.fimg = ttk.Frame(self.fviewSingle)
        self.figure = None
        self.toolbar = None
        self.fimg.pack(side=TOP, expand=True, fill=BOTH)
        self.fviewSingle.pack(side=TOP, fill=BOTH)
        self.fsingleOpt.pack(side=TOP, fill=BOTH)

        self.fmultiOpt=ttk.Frame(self.fcontOptions)
        self.fmoopt=ttk.Frame(self.fmultiOpt)
        self.ve = IntVar()
        ttk.Radiobutton(self.fmoopt, text="All sequences", variable=self.ve, value=0).pack(side=LEFT)
        ttk.Radiobutton(self.fmoopt, text="Any TF", variable=self.ve, value=1).pack(side=LEFT)
        ttk.Radiobutton(self.fmoopt, text="All TF", variable=self.ve, value=2).pack(side=LEFT)
        self.ve.set(0)
        self.fmoopt.pack(side=TOP, fill=X)
        self.fviewMulti = ttk.Frame(self.fmultiOpt)
        self.fviewMulti.pack(side=LEFT, fill=BOTH, expand=True)
        self.fmultiOpt.pack(side=TOP, fill=BOTH)

        self.fcontOptions.add(self.fsingleOpt, text='Single sequence')
        self.fcontOptions.add(self.fmultiOpt, text='Multiple sequences')

        self.fcontOptions.pack(side=TOP, fill=BOTH)

        self.dbbox.bind("<<ComboboxSelected>>", self.updateBox)
        self.confidence.bind("<<ComboboxSelected>>", self.updateAll)
        self.regulator.bind("<<ComboboxSelected>>", self.searchGenes)
        self.fcSequence.bind("<Enter>", self._on_frame_focus)
        self.fcSequence.bind("<Leave>", self._on_frame_lost_focus)
        self.fviewMulti.bind("<Enter>", self._on_frame_focus)
        self.fviewMulti.bind("<Leave>", self._on_frame_lost_focus)

        self.sactree.bind("<ButtonRelease-1>", self._checkItem)
        self.seqtree.bind("<ButtonRelease-1>", self._checkItem)

        self.seq = None
        self.name = None
        self.genome = None
        self.genomeF = None
        self.genes = None
        self.genesF = None
        self.genesTemp = None
        self.genesTempF = None
        self.regulators = None
        self.fileNameF = None
        self.pms = None
        self.pmTargets = None
        self.pmNonTargets = None
        self.reset()
        self.changeInput()
        self.updateBox()
        self.master.mainloop() # após todas especificações da janela

    def _checkItem(self, event):
        tree = self.sactree if(self.fcontainer.index(self.fcontainer.select())==0) else self.seqtree
        if(tree.identify_column(event.x)=="#1"):
            #unchecked b'\xe2\x98\x90' ☐
            #checked b'\xe2\x98\x91' ☑
            #partial b'\xe2\x98\x92' ☒
            genes, genesTemp, search= [self.genes,self.genesTemp,self.sacSearch] if(self.fcontainer.index(self.fcontainer.select())==0) else [self.genesF,self.genesTempF,self.seqSearch]

            if(genes is not None and genesTemp is not None):
                if(tree.identify_region(event.x, event.y)=='heading'):
                    status = tree.heading('check')['text'].encode("utf-8")
                    newStatus = "☑" if(status==b'\xe2\x98\x90' or status==b'\xe2\x98\x92') else "☐"
                    #Mudar todos checks
                    if(search.get().strip()=="" and (self.fcontainer.index(self.fcontainer.select())==1 or self.regulator.get()=='None')): #global
                        genes.loc[:, 'check']=True if(newStatus=="☑") else False
                        genesTemp.loc[:, 'check']=True if(newStatus=="☑") else False
                        for row in tree.get_children():
                            values = tree.item(row)['values']
                            values[0] = newStatus
                            tree.item(row, values=values)
                    else: #local
                        genesTemp.loc[:, 'check']=True if(newStatus=="☑") else False
                        for row in tree.get_children():
                            values = tree.item(row)['values']
                            values[0] = newStatus
                            tree.item(row, values=values)
                            genes.loc[genes['gene']==str(values[1]), 'check']=True if(newStatus=="☑") else False
                else:
                    row = tree.identify_row(event.y)
                    values = tree.item(row)['values']
                    values[0] = "☑" if(tree.item(row)['values'][0].encode("utf-8")==b'\xe2\x98\x90') else "☐"
                    tree.item(row, values=values)

                    genes.loc[genes['gene']==str(values[1]), 'check']=True#if(values[0]=="☑") else False
                    genesTemp.loc[genesTemp['gene']==str(values[1]), 'check']=True if(values[0]=="☑") else False
                    total = len(genesTemp.loc[genesTemp['check']==True])

                    newStatus = "☑" if(total == len(genesTemp)) else ("☐" if total==0 else "☒")
                #Mudar check heading
                tree.heading('check', text=newStatus)

    def _on_frame_focus(self, event):
        if self.OS == "Linux" :
            event.widget.bind_all('<4>', self._on_mousewheel, add="+")
            event.widget.bind_all('<5>', self._on_mousewheel, add="+")
        else: # Windows and MacOS
            event.widget.bind_all("<MouseWheel>", self._on_mousewheel)

    def _on_frame_lost_focus(self, event):
        if self.OS == "Linux" :
            event.widget.unbind_all('<4>')
            event.widget.unbind_all('<5>')
        else: # Windows and MacOS
            event.widget.unbind_all("<MouseWheel>")

    def _on_mousewheel(self, event):
        try:
            delta=0
            if self.OS == 'Linux':
                delta = 2*event.num-9
            elif self.OS == 'Windows':
                delta = (-1)*int((event.delta/120))
            elif self.OS == 'Darwin':
                delta = event.delta

            if(self.fcontOptions.index(self.fcontOptions.select())==0):
                event.widget.xview_scroll(delta, "units")
            else:
                event.widget.yview_scroll(delta, "units")
        except:
            pass

    def Score(self, event=None):
        # create new elements
        targets=[]
        for item in self.vtargets.get():
            aux = item[item.find("(")+1:item.find(")")]
            targets.append(int(aux) if aux!="None" else item.split("/")[0])

        self.targets={}
        if(targets!=[]):
            self.targets = {i:None for i in targets}

        for item in self.vtargets.get():
            x = item.split(" ")
            key = x[1][1:-1]
            key = int(key) if key!="None" else x[0].split("/")[0]
            if(self.targets[key]==None):
                self.targets[key] = (0 if x[2]=="-" else 1)
            else:
                self.targets[key] = 2

        self.maximize = list(self.targets.values())
        self.pmTargets = dbu.getAll(list(self.targets.keys()))[:len(self.targets)]
        targetNames = ["{}/{} ({})".format(pm.StandardName, pm.SystematicName, str(pm.MotifID)) for pm in self.pmTargets]
        self.pmNonTargets = [pm for pm in self.pms if ("{}/{} ({})".format(pm.StandardName, pm.SystematicName, str(pm.MotifID)) not in targetNames and (pm.EC==self.confidence.get() or self.confidence.get()=="All"))]
        if(self.fcontOptions.index(self.fcontOptions.select())==0):
            self.SingleScore()
        elif(self.fcontOptions.index(self.fcontOptions.select())==1):
            self.MultiScore()

    def MultiScore(self,event=None):
        for widget in  self.fviewMulti.winfo_children():
            widget.destroy()

        self.updateGenes()
        self.searchGenes()

        genesTemp = self.genes.loc[self.genes['check']==True] if(self.fcontainer.index(self.fcontainer.select())==0) else self.genesF.loc[self.genesF['check']==True]

        if(len(genesTemp)>0):
            #chamar o canvas
            square=self.master.winfo_height()-(self.fviewMulti.winfo_rooty()-self.master.winfo_rooty())
            cp.canvasPosition(self.fviewMulti, 20, self.fmultiOpt.winfo_width(), 10, seqs=genesTemp,
                targets=[self.targets,self.pmTargets], threshold=float(self.threshold.get()), exist=self.ve.get(), normalized=self.normalizedState.get(), height=square)
        else:
            self.dialog("Plese, select a set.")

    def flatten(self, L):
        for item in L:
            if(isinstance(item,list)):
                yield from self.flatten(item)
            else:
                yield item

    def SingleScore(self, event=None):
        # remove old widgets
        if self.figure:
            self.figure.destroy()
        if self.toolbar:
            self.toolbar.destroy()
        for widget in  self.fcSequence.winfo_children():
            widget.destroy()

        self.name = None
        self.seq = None

        try:
            if(self.fcontainer.index(self.fcontainer.select())==1):
                self.seq = self.sequence.get(1.0,END).rstrip() #rstrip remove quebras de linha #usar RE para verificar integridade da sequência
            elif(self.fcontainer.index(self.fcontainer.select())==2  and len(self.seqtree.get_children())>0):
                #self.seq = self.seqtree.item(self.seqtree.selection())['values'][1]
                self.name = str(self.seqtree.item(self.seqtree.selection())['values'][1])
                self.updateGenes()
                self.searchGenes()
                self.seq = list(self.genesF.loc[self.genesF['gene']==self.name,'sequence'])[0]
            elif(self.fcontainer.index(self.fcontainer.select())==0  and len(self.sactree.get_children())>0):
                self.name = str(self.sactree.item(self.sactree.selection())['values'][1])
                self.updateGenes()
                self.searchGenes()
                self.seq = list(self.genes.loc[self.genes['gene']==self.name,'sequence'])[0]
        except:
            self.dialog("Please, select a sequence.")

        upstream,downstream=[int(self.supstream.get().strip()), -1 if(self.sfgState.get()) else int(self.sdownstream.get().strip())] if(self.fcontainer.index(self.fcontainer.select())==0) else \
            ([int(self.upstream.get().strip()), -1 if(self.fgState.get()) else int(self.downstream.get().strip())] if(self.fcontainer.index(self.fcontainer.select())==2) else [len(self.seq),0])
        if(self.seq is not None):
            cs.canvasSequence(self.seq, self.fcSequence, width=15, height=9, targets=[self.targets,self.pmTargets], outOfRange=self.foorState.get(),
                normalized=self.normalizedState.get(), threshold=float(self.threshold.get()),upstream=upstream,downstream=downstream,margin=60)
            paralogs = dbu.getParalogs(list(self.targets.keys()), type="MotifID", EC=self.confidence.get(), Database=self.dbbox.get()) if(self.paralogState.get()) else []
            paralogs = list(self.flatten(paralogs))

            self.master.update()
            square=min(self.master.winfo_height()-(self.fimg.winfo_rooty()-self.master.winfo_rooty()), self.master.winfo_width()-(self.fimg.winfo_rootx()-self.master.winfo_rootx()))

            plt = Utils.PlotAllScores(self.seq, outOfRange=self.foorState.get(), tfs=self.pmTargets+self.pmNonTargets,
                principalOnly=not(self.allLablesState.get()), save=False, normalized=self.normalizedState.get(),
                paralogs=paralogs, maximize=self.maximize, name=self.name, dpi=self.master.winfo_fpixels('1i'),
                square=square)
            canvas = FigureCanvasTkAgg(plt, self.fimg)
            #self.toolbar = NavigationToolbar2Tk(canvas, self.fimg)
            self.figure = canvas.get_tk_widget()
            self.figure.pack(fill=None)

    def dialog(self,msg):
        messagebox.showinfo("Alerta!" , msg)

    def reset(self):
        self.ltargets.delete(0,END)
        self.lparalogs.config(text='Paralogs = {}')

    def updateAll(self, event=None):
        targets = []
        for item in self.vtargets.get():
            aux = item[item.find("(")+1:item.find(")")]
            targets.append(int(aux) if aux!="None" else item.split("/")[0])

        self.updateSBAll()

        self.lparalogs.config(text='')
        paralogs = [x for x in dbu.getParalogs(targets, type="Order",EC=self.confidence.get(),Database=self.dbbox.get()) if x!=None]
        aux = []
        if(paralogs!=[]):
            if(any(isinstance(i, list) for i in paralogs)):
                paralogs = [item for sublist in paralogs for item in sublist]
            for item in set(paralogs):
                if item not in targets:
                    aux.append(item)
        self.lparalogs.config(text="Paralogs: {"+', '.join(sorted(aux))+"}")

    def updateSBAll(self, event=None):
        self.lAll.delete(0,END)
        search = self.search.get().strip().upper()
        for item in sorted(["{}/{} ({})".format(pm.StandardName,pm.SystematicName,pm.MotifID) for pm in self.pms if pm.EC==self.confidence.get()
            or self.confidence.get()=="All"  if (search in pm.StandardName.upper() or search in pm.SystematicName.upper() or search in str(pm.MotifID) or search=="")]):
            self.lAll.insert(END, item)

    def addTF(self, orientation):
        try:
            item = self.lAll.get(self.lAll.curselection())

            x = item.split(" ")
            exist = False
            targets = []
            for target in self.vtargets.get():
                z = target.split(" ")
                if(x[0]==z[0] and x[1]==z[1] and orientation==z[2]):
                    exist=True

            if(not exist):
                aux = item[item.find("(")+1:item.find(")")]
                targets.append(int(aux) if aux!="None" else item.split("/")[0])
                if(orientation=="+"):
                    self.ltargets.insert(END, item+" +")
                elif(orientation=="-"):
                    self.ltargets.insert(END, item+" -")

            self.lparalogs.config(text='')
            paralogs = [x for x in dbu.getParalogs(targets, type="Order",EC=self.confidence.get(), Database=self.dbbox.get()) if x!=None]
            aux = []
            if(paralogs!=[]):
                if(any(isinstance(i, list) for i in paralogs)):
                    paralogs = [item for sublist in paralogs for item in sublist]
                for item in set(paralogs):
                    if item not in targets:
                        aux.append(item)
            self.lparalogs.config(text="Paralogs: {"+', '.join(sorted(aux))+"}")
        except:
            pass

    def removeTarget(self):
        try:
            selection = self.ltargets.curselection()
            self.ltargets.delete(selection[0])

            targets = []
            for item in self.vtargets.get():
                target = int(item[item.find("(")+1:item.find(")")])
                targets.append(target)

            self.lparalogs.config(text='')
            paralogs = [x for x in dbu.getParalogs(targets, type="Order",EC=self.confidence.get(), Database=self.dbbox.get()) if x!=None]
            aux = []
            if(paralogs!=[]):
                if(any(isinstance(i, list) for i in paralogs)):
                    paralogs = [item for sublist in paralogs for item in sublist]
                for item in set(paralogs):
                    if item not in targets:
                        aux.append(item)
            self.lparalogs.config(text="Paralogs: {"+', '.join(sorted(aux))+"}")
        except:
            pass

    def changeInput(self):
        if(self.fcontainer.index(self.fcontainer.select())==0):
            self.fsac.pack(side=TOP, fill=X)

            if(self.genome==None):
                #carregar reguladores + None
                self.regulators= [None]+dbu.getRegulators()
                self.regulator['values'] = self.regulators
                self.regulator.current(self.regulators.index(None))
                #carregar fasta
                genome = Utils.getGenome("Genome/GCF_000146045.2_R64_genomic.fna")
                #carregar gff
                gff = Utils.getAnnotation("Genome/GCF_000146045.2_R64_genomic.gff")
                self.genome=[genome,gff]

            self.updateGenes()
            self.searchGenes()

    def loadGFF(self):
        fileTypes = ("GFF","*.gff")
        fileTypes=(fileTypes,("all files","*.*"))
        fileNameG=filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes=fileTypes)
        if(fileNameG):
            self.filePathGFF.config(text=fileNameG.split("/")[-1])
        genome = Utils.getGenome(self.fileNameF)
        print("genoma carregado")
        #carregar gff
        gff = Utils.getAnnotation(fileNameG)

        print("gff carregado")
        self.genomeF=[genome,gff]

        self.updateGenes()
        self.searchGenes()

        self.fgOptions.pack(side=TOP, fill=X, expand=YES, padx=2)

    def loadFile(self):
        fileTypes = ("FASTA","*.fasta;*.fna;*.ffn")

        fileTypes=(fileTypes,("all files","*.*"))

        fileName =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes=fileTypes)

        if(fileName):
            self.filePath.config(text=fileName.split("/")[-1])

            if(self.fcontainer.index(self.fcontainer.select())==1):
                if(self.sequence.get(1.0,END).rstrip()!=""):
                    self.sequence.delete(0,END)
                with open(fileName) as f:
                    content = f.readlines()
                self.sequence.insert(END,content)
            elif(self.fcontainer.index(self.fcontainer.select())==2):
                for row in self.seqtree.get_children():
                    self.seqtree.delete(row)
                genes=[]
                with open(fileName, "rU") as handle:
                    i=0
                    for record in SeqIO.parse(handle, "fasta"):
                        genes.append([record.description, str(record.seq), len(record.seq), 0, False])
                        self.seqtree.insert("", i, iid=str(i), values=("☐", record.description,(str(record.seq)[:30]+"..." if len(str(record.seq))>33 else record.seq)))
                        i+=1
                    self.seqtree.selection_set("0")
                self.genesF = pd.DataFrame(genes, columns=["gene", "sequence", "upstream", "downstream", "check"])
                self.fileNameF = fileName
                self.bload.pack(side=LEFT, padx=2)
                self.searchGenes()

    def updateBox(self, event=None):
        self.DBs = dbu.getDBs()
        self.dbbox['values'] = self.DBs['Database'].unique().tolist()
        self.confidence['values'] = ["All" if x is None else x for x in self.DBs['ExpertConfidence'].loc[self.DBs['Database']==self.dbbox.get()].tolist()]
        try:
            self.confidence.current(self.DBs['ExpertConfidence'].loc[self.DBs['Database']==self.dbbox.get()].tolist().index("High"))
        except:
            self.confidence.current(self.DBs['ExpertConfidence'].loc[self.DBs['Database']==self.dbbox.get()].tolist().index(None))
        self.lAll.delete(0,END)
        self.pms = dbu.getAll(Database=self.dbbox.get())
        self.pmTargets = []
        for item in sorted(["{}/{} ({})".format(pm.StandardName,pm.SystematicName,pm.MotifID) for pm in self.pms if pm.EC==self.confidence.get() or self.confidence.get()=="All"]):
            self.lAll.insert(END, item)

    def createDB(self):
        dbname = self.newdb.get().strip()

        if(dbname!="" and dbname!=None):
            if(dbname not in dbu.getDBs()):
                folder = filedialog.askdirectory()
                if(folder!="" and folder!=None):
                    names = dbc.files2DB('TF.db', dbc.SearchFiles(folder,dbname), New=False)
                    #criar as outras tabelas
                    dbc.EmptyInfo('TF.db', names)
                    self.updateBox()
                    self.dialog("Database created!")
            else:
                self.dialog("Nome já existente.")
        else:
            self.dialog("Nome vazio.")

    def fullGene(self, event=None):
        if(self.fcontainer.index(self.fcontainer.select())==0):
            self.sdownstream.configure(state='disabled' if(self.sfgState.get()) else 'normal')
            if(self.sfgState.get()):
                self.sdownstream.delete(0,END)
        elif(self.fcontainer.index(self.fcontainer.select())==2):
            self.downstream.configure(state='disabled' if(self.fgState.get()) else 'normal')
            if(self.fgState.get()):
                self.downstream.delete(0,END)

    def updateGenes(self, event=None):
        if(self.fcontainer.index(self.fcontainer.select())==2 and self.genomeF!=None):
            upstream = int(self.upstream.get().strip())
            downstream = int(self.downstream.get().strip())
            check = self.genesF.iloc[:,4]
            self.genesF = Utils.getSequences(self.genomeF[0], self.genomeF[1], upstream, -1 if(self.sfgState.get()) else downstream)
            self.genesF['check'] = False if(len(check)!=len(self.genesF)) else check

        #pegar genes do FASTA
        if(self.fcontainer.index(self.fcontainer.select())==0):
            upstream = int(self.supstream.get().strip())
            downstream = int(self.sdownstream.get().strip())
            if self.genes is None:
                self.genes = Utils.getSequences(self.genome[0], self.genome[1], upstream, -1 if(self.sfgState.get()) else downstream)
                self.genes['check'] = False
            else:
                check = self.genes.iloc[:,4]
                self.genes = Utils.getSequences(self.genome[0], self.genome[1], upstream, -1 if(self.sfgState.get()) else downstream)
                self.genes['check'] = check


    def searchGenes(self, event=None):
        if(self.fcontainer.index(self.fcontainer.select())==2 and self.genesF is not None):
            search = self.seqSearch.get().strip().upper()
            for row in self.seqtree.get_children():
                self.seqtree.delete(row)
            k=0
            genesTemp=[]
            for i in range(len(self.genesF)):
                name = self.genesF.iloc[i,0]
                sequence = self.genesF.iloc[i,1]
                up = self.genesF.iloc[i,2]
                down = self.genesF.iloc[i,3]
                check = self.genesF.iloc[i,4]
                if(search in name.upper() or search==""):
                    self.seqtree.insert("", i, iid=str(k), values=(('☑' if(check) else '☐'), name,(sequence[:30]+"..." if len(sequence)>33 else sequence)))
                    genesTemp.append([name,sequence,up,down,check])
                    k+=1
            self.genesTempF = pd.DataFrame(genesTemp, columns=["gene", "sequence", "upstream", "downstream", "check"])
            if(k>=1):
                self.seqtree.selection_set("0")
            self.lgenesF['text'] = str(k)+" genes found"
            total = len(self.genesTempF.loc[self.genesTempF['check']==True])
            newStatus = "☑" if(total == len(self.genesTempF)) else ("☐" if total==0 else "☒")
            self.seqtree.heading('check', text=newStatus)

        if(self.fcontainer.index(self.fcontainer.select())==0):
            search = self.sacSearch.get().strip().upper()
            targets = None if(self.regulator.get()=='None') else dbu.getTargets(self.regulator.get().split("/")[1])
            for row in self.sactree.get_children():
                self.sactree.delete(row)
            k=0
            genesTemp=[]
            for i in range(len(self.genes)):
                name = self.genes.iloc[i,0]
                sequence = self.genes.iloc[i,1]
                up = self.genes.iloc[i,2]
                down = self.genes.iloc[i,3]
                check = self.genes.iloc[i,4]
                if(targets is None or name in targets['TargetSys'].tolist() or name in targets['TargetStd'].tolist()) \
                    and (search in name.upper() or search==""): # or search in sequence):
                        #self.sactree.insert("", i, iid=str(k), values=(name,(sequence[:30]+"..."+" ☐ ☑ ☒ " if len(sequence)>33 else sequence)))
                        self.sactree.insert("", i, iid=str(k), values=(('☑' if(check) else '☐'), name,(sequence[:30]+"..." if len(sequence)>33 else sequence)))
                        genesTemp.append([name, sequence,up,down, check])
                        k+=1
            self.genesTemp = pd.DataFrame(genesTemp, columns=["gene", "sequence", "upstream", "downstream", "check"])
            if(k>=1):
                self.sactree.selection_set("0")
            self.lgenes['text'] = str(k)+" genes found"
            total = len(self.genesTemp.loc[self.genesTemp['check']==True])
            newStatus = "☑" if(total == len(self.genesTemp)) else ("☐" if total==0 else "☒")
            self.sactree.heading('check', text=newStatus)

    def saveFile(self):
        self.updateGenes()
        self.searchGenes()
        types = (("FASTA files","*.faa;*.frn;*.ffn;*.fna;*.fasta"),("All files","*.*"))
        f = filedialog.asksaveasfile(mode='w', defaultextension=".fasta", filetypes=types)
        if f is None:
            return
        if(self.fcontainer.index(self.fcontainer.select())==2):
            for i in range(len(self.genesTempF)):
                f.write(">"+self.genesTempF.iloc[i,0]+"\n")
                f.write(self.genesTempF.iloc[i,1]+"\n")
            f.close()
        if(self.fcontainer.index(self.fcontainer.select())==0):
            for i in range(len(self.genesTemp)):
                f.write(">"+self.genesTemp.iloc[i,0]+"\n")
                f.write(self.genesTemp.iloc[i,1]+"\n")
            f.close()

    def validate(self, P, s, S, W, V, dv, type):
        if(type=="Int"):
            if(P.strip()=="" and V=="focusout"):
                self.master.nametowidget(W).delete(0,END)
                self.master.nametowidget(W).insert(0,dv)
                return True
            try:
                if(P!=""):
                    x = int(P)
                return True
            except:
                return False
        elif(type=="Float"):
            if(P.strip()=="" and V=="focusout"):
                self.master.nametowidget(W).delete(0,END)
                self.master.nametowidget(W).insert(0,dv)
                return True
            try:
                if(P!=""):
                    x = float(P)
                return True
            except:
                return False
        elif(type=="DNA"):
            return (S in 'actgnxACTGNX')

if __name__ == '__main__':
    if(not os.path.isfile("TF.db")):
        dbc.main()
    InterfaceView()
