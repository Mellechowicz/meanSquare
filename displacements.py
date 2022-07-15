#!/usr/bin/python3

import numpy as np
import yaml

import argv

def symmetric_matix(vector):
    """
    FROM: Phonopy 6 vector written as [xx, yy, zz, yz, xz, xy]
    TO:   numPy 3x3 matrix
    """
    return np.matrix([[vector[0], vector[5], vector[4]], 
                      [vector[5], vector[1], vector[3]], 
                      [vector[4], vector[3], vector[2]]])


class YAMLFiles:
    """
    Holding opened/closed YAML files
    """
    def __init__(self,*args,**kwargs):
        self.names = args
        self.files = []

        if 'dry_run' in kwargs:
            if kwargs['dry_run']:
                return
        self.open()

    def open(self):
        if self.files:
            raise RuntimeError('Files already opened!',*self.names)
        for name in self.names:
            try:
                yamlfile = open(name)
            except FileNotFoundError as exc:
                print(exc,end='')
                print(' ................omiting')
                yamlfile.close()
                continue
            try:
                self.files.append(yaml.safe_load(yamlfile))
            except yaml.YAMLError as exc:
                print(exc,end='')
                print(' ................omiting')
            yamlfile.close()

    def close(self):
        for idx,file in reversed(list(enumerate(self.files))):
            del self.files[idx]
