#!/usr/bin/python3

import numpy as np
import yaml

import argv
from displacements import symmetric_matix,YAMLFiles
from auxplot import NicePlot

if __name__ != '__main__':
    exit(255)

cmdLine = argv.Options()

if cmdLine('qha'):
    allDisplacements = []
    fileVolumes = np.loadtxt(cmdLine('e_v'))[:,0]
    volumeTemp  = np.loadtxt(cmdLine('volume_temperature'))

files_read = YAMLFiles(*cmdLine('files'))

for idx,databaseConfig in enumerate(files_read.files):
    yaml_file_name = files_read.names[idx][:-5]
    atoms = cmdLine('indexes')
    NATOM   = databaseConfig['natom']
    TEMPS   = []
    mean_sq = [[] for a in atoms]
    for T in databaseConfig['thermal_displacement_matrices']:
        TEMPS.append(T['temperature'])
        for i,atom in enumerate(atoms):
            mean_sq[i].append(
                        np.sum(
                  np.linalg.eig(
                 symmetric_matix(
                                 T['displacement_matrices'][atom]
                                )
                               )[0] #the eigenvalues <=> displacements along the principal axes
                              )
                             )
    displacements = np.transpose(np.array([TEMPS, *mean_sq]))

    np.savetxt(yaml_file_name+'.txt',displacements)

    if cmdLine('plot'):
        myplot = NicePlot(uselatex=cmdLine('latex'),atoms=atoms,labels=cmdLine('labels'))
        myplot.plot(displacements)

        if not cmdLine('qha'):
            myplot.set_title()
        else:
            myplot.set_title(fileVolumes[idx])

        myplot.savefig(yaml_file_name+'.png',dpi=cmdLine('dpi'))

    if cmdLine('qha'):
        allDisplacements.append(displacements)

if cmdLine('qha'):
    allDisplacements = np.array(allDisplacements)
    if fileVolumes.shape[0] != allDisplacements.shape[0]:
        raise IndexError("The data in "+cmdLine('e_v')+" does not match provided files!")
    from scipy.interpolate import interp1d
    temperatures     = allDisplacements[0,:,0]
    qhaDisplacements = []
    for i,T in enumerate(temperatures):
        volume        = volumeTemp[volumeTemp[:, 0] == T][0,1]
        displacements = allDisplacements[:,i,1:]
        interpolation = np.array([ interp1d(fileVolumes,displacements[:,j])(volume) for j,_ in enumerate(displacements[0,:]) ])
        qhaDisplacements.append([T,*interpolation,volume])

    qhaDisplacements = np.array(qhaDisplacements)
    np.savetxt("qha_displacements.txt",qhaDisplacements,header="# temperature volume  "+''.join([label+"  " for label in myplot.labels]))

    if cmdLine('plot'):
        myplot = NicePlot(uselatex=cmdLine('latex'),atoms=atoms,labels=cmdLine('labels'))
        myplot.set_title(r"Mean-square displacement (Quasiharmonic approximation)")
        myplot.plot(qhaDisplacements[:,:-1])
        myplot.twinx(qhaDisplacements[:,0],qhaDisplacements[:,-1])
        myplot.savefig('qha_displacement.png',dpi=cmdLine('dpi'))
