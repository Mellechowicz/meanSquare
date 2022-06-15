# -*- coding: utf-8 -*-
import argparse as ap
import numpy as np
import re
from sys import argv

class Error:
    access = 123

class Options:
    def __init__(self):
        self.parser = ap.ArgumentParser(description='Extracting mean-square displacement along the principal axes')
        self.add_arguments_input()
        self.add_arguments_output()
        self.add_arguments_matplotlib()

        self.opt = self.parser.parse_args(argv[1:])
        self.prepare()

    def add_arguments_input(self):
        self.parser.add_argument('--files','-f',default=["thermal_displacement_matrices.yaml"],nargs='+',
                                 help="a set of thermal_displacement_matrices.yaml files")
        self.parser.add_argument('--indexes','-i',default=["-1"],nargs='+',
                                 help="index of atom in supercell (staring from 0), default: last atom")
        self.parser.add_argument('--e-v','-V',default="e-v.dat",nargs=1,
                                 help="phonopy \'e-v.dat\' file corresponding to thermal displacement files")
        self.parser.add_argument('--volume-temperature','-v',default="volume-temperature.dat",nargs=1,
                                 help="phonopy \'volume-temperature.dat\' file")

    def add_arguments_output(self):
        self.parser.add_argument('--plot','-p', action='store_true',
             help='creates a plot of displacement vs. temperature')
        self.parser.add_argument('--qha','-q', action='store_true',
             help='use phonPy output of quasiharmonic approximation ')

    def add_arguments_matplotlib(self):
        self.parser.add_argument('--labels','-l',default=[],nargs='+',
                                 help="a set of labels in the resultant figure")
    def prepare(self):
        self.fix_files()
        self.fix_indexes()

    def fix_files(self):
        self.files = []
        for f in self.opt.__dict__['files']:
            if f[-4:] == 'yaml':
                self.files.append(f)
            else:
                print("Ignoring file %s"%f)
        self.opt.__dict__['files'] = self.files

    def fix_indexes(self):
        self.indexes= []
        for i in self.opt.__dict__['indexes']:
            try:
                self.indexes.append(int(i))
            except ValueError:
                print("Ignoring index %s"%i)
        self.opt.__dict__['indexes'] = self.indexes

    def __call__(self, key):
        return self.opt.__dict__[key]
