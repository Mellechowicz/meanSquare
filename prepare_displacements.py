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

cmdLine = argv.Options()

print(cmdLine)
print(cmdLine('dpi'))
if cmdLine('qha'):
    allDisplacements = []
    fileVolumes = np.loadtxt(cmdLine('e_v'))[:,0]
    volumeTemp  = np.loadtxt(cmdLine('volume_temperature'))

for idx,yaml_file in enumerate(cmdLine('files')):
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

            if cmdLine('plot'):
                import matplotlib.pyplot as plt
                if cmdLine('latex'):
                    plt.rcParams.update({ "text.usetex": True,
                                          "font.family": "serif"})
                fig,ax = plt.subplots()
                labels = cmdLine('labels')
                diff   = len(atoms)-len(labels) 
                if diff > 0:
                    for i in range(len(atoms)-diff,len(atoms)):
                        labels.append("atom #%d"%atoms[i])
                labels = labels[:len(atoms)]
                for i,label in enumerate(labels):
                    ax.plot(displacements[:,0],displacements[:,i+1],label=label)
                if not cmdLine('qha'):
                    ax.set_title("Mean-square displacement (constant volume)")
                else:
                    if cmdLine('latex'):
                        ax.set_title(r"Mean-square displacement at $V$= "+"%9.0f"%fileVolumes[idx]+r" (\AA{}$^3$)")
                    else:
                        ax.set_title("Mean-square displacement (V=%9.0f A^3)"%fileVolumes[idx])
                if cmdLine('latex'):
                    ax.set_xlabel(r"temperature, $T$ (K)")
                    ax.set_ylabel(r"displacement, $d^2$ (\AA{}$^2$)")
                else:
                    ax.set_xlabel("temperature, T (K)")
                    ax.set_ylabel("displacement, d^2 (A^2)")
                ax.grid()
                ax.legend()
                plt.savefig(yaml_file[:-5]+'.png',dpi=cmdLine('dpi'))

            if cmdLine('qha'):
                allDisplacements.append(displacements)
        except yaml.YAMLError as exc:
            print(exc)

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
        qhaDisplacements.append([T,volume,*interpolation])
    qhaDisplacements = np.array(qhaDisplacements)
    if cmdLine('plot'):
        import matplotlib.pyplot as plt
        if cmdLine('latex'):
            plt.rcParams.update({ "text.usetex": True,
                                  "font.family": "serif"})
        fig,ax = plt.subplots()
        for i,label in enumerate(labels):
            ax.plot(qhaDisplacements[:,0],qhaDisplacements[:,i+2],label=label)
        if cmdLine('latex'):
            ax.set_title (r"Mean-square displacement (Quasiharmonic approximation)")
            ax.set_xlabel(r"temperature, $T$ (K)")
            ax.set_ylabel(r"displacement, $d^2$ (\AA{}$^2$)")
        else:
            ax.set_title( "Mean-square displacement (Quasiharmonic approximation)")
            ax.set_xlabel("temperature, T (K)")
            ax.set_ylabel("displacement, d^2 (A^2)")

        ax2 = ax.twinx()
        ax2.plot(qhaDisplacements[:,0],qhaDisplacements[:,1],'--',c='k',label='volume')
        if cmdLine('latex'):
            ax2.set_ylabel(r"volume, $V$ (\AA{}$^3$)")
        else:
            ax2.set_ylabel("volume, V (A^3)")
        ax.plot([],[],'--',c='k',label='cell volume') # Just for the legend
        limits = ax2.get_ylim()
        limits = (1.3*limits[0]-0.3*limits[1],limits[1]) # so that volume doesn't overlap with the displacements
        ax2.set_ylim(limits)

        ax.grid()
        ax.legend()
        plt.savefig('qha_displacement.png',dpi=cmdLine('dpi'))
