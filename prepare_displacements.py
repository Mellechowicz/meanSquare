#!/usr/bin/python3

import numpy as np
import yaml

import argv
# argv[1]: mean_sqr_disp yaml file
# argv[2]: atom of interest

def symmetric_matix(vector):
    """
    FROM: Phonopy 6 vector written as [xx, yy, zz, yz, xz, xy]
    TO:   numPy 3x3 matrix
    """
    return np.matrix([[vector[0], vector[5], vector[4]], 
                      [vector[5], vector[1], vector[3]], 
                      [vector[4], vector[3], vector[2]]])

cmdLine = argv.Options()

for yaml_file in cmdLine('files'):
    with open(yaml_file) as file:
        atoms = cmdLine('indexes')
        try:
            databaseConfig = yaml.safe_load(file)   
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

            np.savetxt(yaml_file[:-5]+'.txt',displacements)

            if not cmdLine('plot'):
                break
    
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots()
            labels = cmdLine('labels')
            diff   = len(atoms)-len(labels) 
            if diff > 0:
                for i in range(len(atoms)-diff,len(atoms)):
                    labels.append("atom #%d"%atoms[i])
            labels = labels[:len(atoms)]
            for i,label in enumerate(labels):
                ax.plot(displacements[:,0],displacements[:,i+1],label=label)
                ax.set_title("Mean-square displacement")
                ax.set_xlabel("temperature, T (K)")
                ax.set_ylabel("displacement, d^2 (A^2)")
            ax.grid()
            ax.legend()
            plt.savefig(yaml_file[:-5]+'.png',dpi=300)
        except yaml.YAMLError as exc:
            print(exc)
