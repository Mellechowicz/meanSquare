#!/usr/bin/python3

import matplotlib.pyplot as plt
import importlib
from numbers import Number

class NicePlot:
    temperature = {
            False: "temperature, T (K)",
            True: r"temperature, $T$ (K)"}
    displacement = {
            False: "displacement, d^2 (A^2)",
            True: r"displacement, $d^2$ (\AA{}$^2$)"}
    titlePrefix = {
            False: "Mean-square displacement at V = ",
            True: r"Mean-square displacement at $V$ = "}
    titleSuffix = {
            False: " (A^3)",
            True: r" (\AA{}$^3$)"}

    def __init__(self,uselatex=False,atoms=[],labels=[]):
        importlib.reload(plt)

        self.uselatex = uselatex
        if self.uselatex:
            plt.rcParams.update({ "text.usetex": True,
                                  "font.family": "serif"})

        self.set_atoms(atoms)
        self.set_labels(labels)
        self.fig, self.ax = plt.subplots()

    def set_atoms(self,atoms):
        self.atoms = atoms

    def set_labels(self,labels):
        self.labels = labels

    def fix_labels(self):
        diff   = len(self.atoms)-len(self.labels) 
        if diff > 0:
            for i in range(len(self.atoms)-diff,len(self.atoms)):
                self.labels.append("atom %d"%self.atoms[i])
        self.labels = self.labels[:len(self.atoms)]

    def single_plot(self,x,y,label):
        self.ax.plot(x,y,label=label)

    def plot(self,data):
        self.fix_labels()
        for i in range(1,data.shape[1]): #iterate over columns
            self.ax.plot(data[:,0],data[:,i],label=self.labels[i-1])

        self.set_axes()
        self.ax.grid()
        self.ax.legend()

    def set_axes(self):
        self.ax.set_xlabel(self.temperature[self.uselatex])
        self.ax.set_ylabel(self.displacement[self.uselatex])

    def set_title(self,volume_maybe=-1.0):
        if isinstance(volume_maybe,str):
            self.ax.set_title(volume_maybe)
            return
        title_string = self.titlePrefix[self.uselatex]
        if isinstance(volume_maybe,Number):
            title_string += "%9.0f"%volume_maybe \
                           +self.titleSuffix[self.uselatex]
        else:
            title_string += "const"
        self.ax.set_title(title_string)

    def savefig(self,name,dpi=300):
        plt.savefig(name,dpi=dpi)

    def twinx(self,x,y,label='volume',y2limits=None):
        ax2 = self.ax.twinx()
        ax2.plot(x,y,'--',c='k',label=label)
        ax2.set_ylabel(r"volume, V" + self.titleSuffix[self.uselatex])
        self.ax.plot([],[],'--',c='k',label='cell volume') # Just for the legend
        if y2limits is None:
            limits = ax2.get_ylim()
            limits = (1.3*limits[0]-0.3*limits[1],limits[1]) # so that volume doesn't overlap with the displacements
        else:
            limits = y2limits
        ax2.set_ylim(limits)
