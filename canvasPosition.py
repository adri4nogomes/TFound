from tkinter import *
import tkinter as tk
from tkinter import font
import random
import pandas as pd
import Utils
import numpy as np
import PositionMatrix as pm
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

COLORS = ['dark slate gray', 'dim gray', 'slate gray',
'gray','midnight blue', 'navy', 'cornflower blue', 'dark slate blue',
'slate blue', 'medium slate blue', 'medium blue', 'royal blue',  'blue',
'dodger blue', 'deep sky blue', 'sky blue', 'steel blue',
'powder blue', 'pale turquoise', 'dark turquoise', 'medium turquoise', 'turquoise',
'cyan', 'cadet blue', 'medium aquamarine', 'aquamarine', 'dark green', 'dark olive green',
'dark sea green', 'sea green', 'medium sea green', 'pale green', 'spring green',
'lawn green', 'medium spring green', 'green yellow', 'lime green', 'yellow green',
'forest green', 'olive drab', 'dark khaki', 'khaki', 'pale goldenrod',
'yellow', 'gold', 'goldenrod', 'dark goldenrod', 'rosy brown',
'indian red', 'saddle brown', 'sandy brown','dark salmon', 'salmon', 'orange', 'dark orange',
'coral', 'tomato', 'orange red', 'red', 'hot pink', 'deep pink', 'pink',
'pale violet red', 'maroon', 'medium violet red', 'violet red',
'medium orchid', 'dark orchid', 'dark violet', 'blue violet', 'purple', 'medium purple',
'thistle', 'SlateBlue1','SlateBlue4', 'RoyalBlue1', 'RoyalBlue4', 'blue4',
'DodgerBlue2', 'DodgerBlue3', 'DodgerBlue4', 'SteelBlue1', 'SteelBlue4', 'DeepSkyBlue2', 'DeepSkyBlue4',
'SkyBlue1', 'SkyBlue4', 'SlateGray1', 'SlateGray4', 'PaleTurquoise1', 'PaleTurquoise4', 'CadetBlue1',
'CadetBlue4', 'turquoise1', 'turquoise4','cyan4', 'DarkSlateGray1', 'DarkSlateGray4',
'aquamarine2', 'aquamarine4', 'DarkSeaGreen1', 'DarkSeaGreen4', 'SeaGreen1', 'SeaGreen3',
'PaleGreen1', 'PaleGreen4', 'SpringGreen2', 'SpringGreen4', 'green4', 'chartreuse2',
'chartreuse4','OliveDrab1', 'OliveDrab4', 'DarkOliveGreen1','DarkOliveGreen4', 'khaki1', 'khaki4',
'yellow4','gold4', 'goldenrod1', 'goldenrod4', 'DarkGoldenrod1', 'DarkGoldenrod2', 'DarkGoldenrod4',
'RosyBrown1', 'RosyBrown4', 'IndianRed1', 'IndianRed4', 'sienna1', 'sienna4', 'burlywood1',
'burlywood4', 'wheat1', 'wheat4', 'tan1','tan4', 'chocolate1', 'firebrick1',  'firebrick4',
'brown1', 'brown4', 'salmon1','salmon4', 'orange4', 'DarkOrange1',
'DarkOrange4', 'coral1', 'coral4', 'tomato2', 'tomato4', 'OrangeRed2','OrangeRed4', 'red4',
'DeepPink4','HotPink1', 'HotPink4', 'pink1', 'pink4', 'PaleVioletRed1',
'PaleVioletRed4', 'maroon1', 'maroon4', 'VioletRed1', 'VioletRed4','magenta2', 'magenta4', 'orchid1',
'orchid4', 'plum1','plum4', 'MediumOrchid1', 'MediumOrchid4', 'DarkOrchid1', 'DarkOrchid4','purple1',
'purple4', 'MediumPurple1', 'MediumPurple4', 'thistle1', 'thistle4']

class canvasPosition():
    def __init__(self, frame, margin, width, rad, seqs, targets, threshold, exist, normalized=True, height=0):
        self.margin = margin
        self.width = width - 16 #-scrollbar
        self.rad = rad
        self.targets = targets[0]
        self.seqs = seqs
        self.tfs = targets[1]
        self.threshold = threshold
        self.exist = exist
        self.normalized=normalized

        self.position = Frame(frame)
        self.canvas = Canvas(self.position, bg='white')
        bar=Scrollbar(self.position, orient=VERTICAL)
        bar.pack(side=RIGHT, fill=Y)
        bar.config(command=self.canvas.yview, width=10)
        self.canvas.config(yscrollcommand=bar.set)
        self.canvas.pack(side=LEFT, expand=True, fill=BOTH)
        self.position.pack(side=TOP,fill=BOTH, expand=True)
        #criar círculos principais
        y, colors = self.circles()

        #criar linhas das Sequências
        x,y,maxUp,maxDown = self.circleLine(y, colors)
        #criar círculos nas linhas das sequências
        #criar régua
        self.scalarLine(x,y,maxUp,maxDown)

        self.canvas.config(width=self.width,height=y+3*self.rad if height==0 else height)
        self.canvas.config(scrollregion=(0,0,0,y+3*self.rad))

    def circles(self):
        global COLORS
        width=self.width
        margin = self.margin
        tfs=self.tfs
        targets=self.targets
        colors = [None]
        if(len(tfs)>=1):
            colors = random.sample(COLORS, len(tfs))
        x=self.margin
        y=self.margin
        if(len(tfs)>0):
            for i in range(len(tfs)):
                x,y = self.circle(x,y,colors[i],list(self.targets.values())[i],tfs[i].StandardName)
                x+=3*self.rad
            return [y+3*self.rad,colors]
        return [y,colors]


    def circle(self, x, y, color, type, title=None, mini=False):
        xf, yf = [x, y]
        r = self.rad if not mini else self.rad/2
        bbox = None
        text_item = None
        if(title!=None):
            while(True):
                text_item = self.canvas.create_text(x+1.5*r,y,text=title, anchor=W, fill="black", font=('Arial', -(r*2), 'bold'))
                bbox = self.canvas.bbox(text_item)
                if(bbox[2]<self.width-self.margin):
                    xf=bbox[2]
                    yf=y
                    break
                else:
                    self.canvas.delete(text_item)
                    x=self.margin
                    y+=self.margin+1*self.rad
                    text_item = None

        if(not mini):
            self.canvas.create_oval(x-(r if (type==2 or title!=None) else (2*r if type==0 else 0)),y-r,x+(r if (type==2 or title!=None) else (2*r if type==1 else 0)),y+r, fill=(color if type==2 else 'white'), width=1)
        if(type!=2):
            self.canvas.create_arc(x-(r if (type==2 or title!=None) else (2*r if type==0 else 0)),y-r,x+(r if (type==2 or title!=None) else (2*r if type==1 else 0)),y+r, start=0 if type==1 else 180, extent=180, fill=color)
        return [xf,yf]

    def circleLine(self, y, colors):
        underline = False
        maxUp, maxDown = max(self.seqs.loc[:,'upstream']), max(self.seqs.loc[:,'downstream'])
        fullLen = maxUp+maxDown
        maxName = max(len(x) for x in self.seqs.loc[:,'gene'])
        x = (maxName+4)*self.rad
        total = self.width-self.margin-x
        nseq = 0

        #fig = mpl.figure.Figure(figsize=(1, 1))
        #ax = fig.add_axes([0, 0, 1, 1])
        concatScoreS = None
        concatScoreA = None
        for seq in range(len(self.seqs)):
            scoresS = [None]*(len(self.tfs))
            scoresA = [None]*(len(self.tfs))
            exist = [True]*(len(self.tfs)*2)
            for i in range(len(self.tfs)):
                if(list(self.targets.values())[i] in [1,2]):
                    scoresS[i] = Utils.getScorePostion(self.seqs.iloc[seq,1],self.tfs[i], self.threshold, self.normalized, outOfRange=False, reversePosition=False)
                    exist[i] = False if scoresS[i].size==0 else True

                if(list(self.targets.values())[i] in [0,2]):
                    scoresA[i] = Utils.getScorePostion(self.seqs.iloc[seq,1],self.tfs[i], self.threshold, self.normalized, outOfRange=False, reversePosition=True)
                    exist[i+len(self.tfs)] = False if scoresA[i].size==0 else True

                exist[i+len(self.tfs)] = exist[i] if(list(self.targets.values())[i]==1) else exist[i+len(self.tfs)]
                exist[i] = exist[i+len(self.tfs)] if(list(self.targets.values())[i]==0) else exist[i]

            if(self.exist==0 or (self.exist==1 and any(exist)) or (self.exist==2 and all(exist))):
                nseq+=1
                f = tk.font.Font(size=-2*self.rad, weight='bold', underline=(1 if(underline) else 0))
                self.canvas.create_text(self.rad,y,text=self.seqs.iloc[seq,0], anchor=W, fill="black", font=f)

                upstream = self.seqs.iloc[seq,2]
                downstream = self.seqs.iloc[seq,3]
                self.canvas.create_line(x+(self.width-self.margin-x)*((maxUp-upstream)/fullLen), y, x+(self.width-self.margin-x)*((maxUp+downstream)/fullLen), y, fill='#7F7F7F', width=2, arrow=tk.LAST, capstyle=tk.ROUND)
                if(maxDown>0):
                    self.transversal(self.width-self.margin-total*(maxDown/fullLen), y)
                for i in range(len(self.tfs)):
                    if(scoresS[i] is not None):
                        for s in scoresS[i]:
                            self.circle(x+total*(int(s[1])/fullLen), y, color=colors[i], type=1, mini=True)
                    if(scoresA[i] is not None):
                        for s in scoresA[i]:
                            self.circle(x+total*(int(s[1])/fullLen), y, color=colors[i], type=0, mini=True)
                y+=3*self.rad
        self.canvas.create_text(self.rad,y,text=str(nseq)+" genes", anchor=W, fill="#7F7F7F", font=('Arial', -(int(self.rad*1.2)), 'bold'))
        return [x,y, maxUp, maxDown]

    def scalarLine(self,x, y, upstream, downstream, interv=5, tam=10):
        self.canvas.create_line(x, y, self.width-self.margin, y, fill='#7F7F7F', width=2, capstyle=tk.ROUND)
        i = 0
        total = upstream+downstream
        length = self.width-self.margin-x
        while(i<=total):
            if(abs(i-upstream)>(total/(4*interv))):
                self.transversal(x+length*(i/total), y, ("+" if i>upstream else "")+str(int(i-upstream)))
            i+=total/interv
        self.transversal(x+length*(upstream/total), y, "ATG")

    def transversal(self, x, y, text=None):
        self.canvas.create_line(x, y+5, x, y-5, fill='#7F7F7F', width=2, capstyle=tk.ROUND)
        if(text!=None):
            self.canvas.create_text(x, y+10, text=text, anchor=N, fill="#7F7F7F", font=('Arial', -(int(self.rad*1.2)), 'bold'))
