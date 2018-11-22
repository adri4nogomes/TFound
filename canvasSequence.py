from tkinter import *
import tkinter as tk
from tkinter import ttk
from tkinter import font
import Utils
import PositionMatrix as pm
import numpy as np
import random

COLORS = ['dark slate gray', 'dim gray', 'slate gray',
'gray','midnight blue', 'navy', 'cornflower blue', 'dark slate blue',
'slate blue', 'royal blue',  'blue',
'dodger blue', 'deep sky blue', 'sky blue', 'steel blue',
'powder blue', 'pale turquoise', 'dark turquoise', 'turquoise',
'cyan', 'cadet blue', 'aquamarine', 'dark green', 'dark olive green',
'dark sea green', 'sea green', 'pale green', 'spring green',
'lawn green', 'green yellow', 'lime green', 'yellow green',
'forest green', 'olive drab', 'dark khaki', 'khaki', 'pale goldenrod',
'yellow', 'gold', 'goldenrod', 'dark goldenrod', 'rosy brown',
'indian red', 'saddle brown', 'sandy brown','dark salmon', 'salmon', 'orange', 'dark orange',
'coral', 'tomato', 'orange red', 'red', 'hot pink', 'deep pink', 'pink',
'pale violet red', 'maroon', 'violet red',
'dark orchid', 'dark violet', 'blue violet', 'purple',
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

class canvasSequence():
    def __init__(self, seq, frame, width=15, height=7, targets=None, threshold=0.8, normalized=True, outOfRange=False, upstream=0, downstream=0, margin=0):
        self.tab = str.maketrans("atgcATGC", "tacgTACG")
        self.width = width
        self.height = height
        self.seq = seq
        self.targets = targets
        self.threshold = threshold
        self.normalized = normalized
        self.outOfRange = outOfRange
        self.upstream = upstream
        self.downstream = downstream
        self.margin = margin

        self.canvas = Canvas(frame, scrollregion=(0,0,len(self.seq)*width+2*margin,0),bg="white")
        bar=Scrollbar(frame, orient=HORIZONTAL)
        bar.pack(side=TOP, fill=X)
        bar.config(command=self.canvas.xview, width=10)
        self.canvas.config(width=500,height=height+5.5*width+5*width/5)
        self.canvas.config(xscrollcommand=bar.set)
        self.canvas.pack(side=LEFT, expand=True, fill=BOTH)

        self.sequence(seq)

        self.arrow()

    def color(self, base):
        colors = {'A':'#FD3333', 'T':'#35C430', 'C': '#3044C4', 'G':'yellow'}
        if(base in colors):
            return colors[base]
        else:
            return "black"

    def arrow(self):
        global COLORS
        tfs=self.targets[1]
        targets=self.targets[0]
        colors = [None]
        if(len(tfs)>1):
            colors = random.sample(COLORS, len(tfs))
        #loop por PWMS
        for i in range(len(list(targets.keys()))):
            #loop por resultados do PWM superiores ao threshold
            if(list(targets.values())[i] in [1,2]):
                scores = Utils.getScorePostion(self.seq, tfs[i], self.threshold, self.normalized, outOfRange=self.outOfRange, reversePosition=False)
                for s in scores:
                    self.line(int(s[1]),len(tfs[i].pm), name=tfs[i].StandardName, score=s[0], sense=True,color=colors[i])

            if(list(targets.values())[i] in [0,2]):
                scores = Utils.getScorePostion(self.seq, tfs[i], self.threshold, self.normalized, outOfRange=self.outOfRange, reversePosition=True)
                for s in scores:
                    self.line(int(s[1]),len(tfs[i].pm), name=tfs[i].StandardName, score=s[0], sense=False,color=colors[i])

    def line(self, bi, n, name=None, score=None, sense=True, arrowSense=True,color='black'):
        height = self.height
        width = self.width
        ini = self.margin+bi*width if sense else self.margin+(bi+1)*width
        end = self.margin+bi*width+n*width if sense else self.margin+(bi+1)*width-n*width
        space = width/5
        lineDist = (height+1.5*width+1.25*space if sense else height+3.5*width+4.5*space)
        arrowshape = [x*(width/20) for x in [10,10,2.5]]
        self.canvas.create_line(ini, lineDist, end, lineDist, width=space, arrow=(tk.LAST if arrowSense else tk.FIRST), arrowshape=arrowshape, fil=color)
        recDist = lineDist+(1.75*space if sense else -(1.5*space+width))
        self.canvas.create_rectangle(ini, recDist, end, recDist+width, width=width/5, outline=color)
        self.canvas.create_text(((ini+end)/2),(height if sense else height+3.5*width+4.5*space),anchor=N,text=(name if name!=None else ""), font=('Helvetica', -width,'bold'))
        self.canvas.create_text(((ini+end)/2),(height+width if sense else height+4.5*width+4.5*space),anchor=N,text="({0:.3f})".format(score), font=('Helvetica', -int(round(width/2)),'bold'))

    def create_box(self, letter, xi, yi, width, wfont, height=None, border=True,color=False):
        xf = xi+width
        yf = yi+width if height==None else yi+height
        self.canvas.create_rectangle(xi, yi, xf, yf ,fill=(self.color(letter) if color else None),width=(1 if border else 0))
        self.canvas.create_text((xi+xf)/2,width/15+(yf+yi)/2,text=letter, fill="black", font=('Helvetica', -wfont, 'bold'))

    def sequence(self, antisense=True):
        height = self.height
        width = self.width
        for j in range(len(self.seq)):
            ini = self.margin+j*width
            #("-"+str(j+1) if((j+1)%10==0 or (j+1)==0) #numeração padrão
            position = ("-"+str(self.upstream-j) if((j)%10==0 or (j+1)==self.upstream) else ".") if j<self.upstream else ("+"+str(j+1-self.upstream) if((j+1)%10==0 or (j-self.upstream)==0) else ".")
            self.create_box(position, ini, 0,width,height, height, border=False)
            self.create_box(self.seq[j], ini, height+1.5*width+3*width/5,width,width)
            if(antisense):
                self.create_box(self.seq[j].translate(self.tab), ini, height+2.5*width+3*width/5,width,width)
