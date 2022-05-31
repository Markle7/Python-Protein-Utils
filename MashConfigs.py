def show_exception_and_exit(exc_type, exc_value, tb):
    import traceback
    traceback.print_exception(exc_type, exc_value, tb)
    input("Press key to exit.")
    sys.exit(-1)


import numpy as np
import os
import os.path as path
import sys

sys.excepthook = show_exception_and_exit

def MashConfigurations(argv):
    print("Mashing all conformations together")
    if (len(argv) != 2):
        raise ValueError("Specified more or less than 2 files. Required 2 files. %.conf and %.pdb files")
        
    #Set docking and pdb file paths
    confPath = [x.replace("\\","/") for x in argv if ".conf" in x][0]
    pdbPath = [x.replace("\\","/") for x in argv if ".pdb" in x][0]

    if (confPath == "" or pdbPath == ""):
        raise ValueError("Missing one file")


    if True:
        #THIS PROGRAM IS TAILORED TO BE USED WITH QUIMERA_KV12_21 PDB.
        #FOR THIS PDB, THE CHAINS ASSOCIATION WITH SHAKER IS AS FOLLOWS: (SHAKER CHAIN ID, CHIMERA CHAIN ID)
        #(A,B) ; (B,F) ; (C,G) ; (D, H)
        #To make this program flexible we define an array that transforms (A,B,C,D) to anyother chain layout

        #Chain layout for Shaker PDB
        chainLayout = ["A","B","C","D"]

        #Chain layout for Quimera PDB
        #chainLayout = ["B", "F", "G", "H"]

        #The following algorithm calculates the rotation matrix required to align chain X onto chain D
        PDBfile = open(pdbPath,"r")
        oPDB = PDBfile.read().split("\n")
        PDBfile.close()
        chains = []
        chains += [[a for a in oPDB if len(a) > 22 and a[0:4] == "ATOM" and a[21] == chainLayout[0]]]
        chains += [[a for a in oPDB if len(a) > 22 and a[0:4] == "ATOM" and a[21] == chainLayout[1]]]
        chains += [[a for a in oPDB if len(a) > 22 and a[0:4] == "ATOM" and a[21] == chainLayout[2]]]
        chains += [[a for a in oPDB if len(a) > 22 and a[0:4] == "ATOM" and a[21] == chainLayout[3]]]
        #R defines a list of the rotation matrices to take chainX to chainD. In order: A, B, C
        R = []
        C = []
        T = []
        
        pdbAllAtoms = []
        for i in range(4):
            for j in range(len(chains[i])):
                pdbAllAtoms += [np.asarray([chains[i][j][30:38],chains[i][j][38:46],chains[i][j][46:54]],dtype="float64")]
        pdbOrigin = np.sum(pdbAllAtoms,axis=0)/len(pdbAllAtoms)
        
        for i in range(3):
            print("\tReading chain " + str(i) + " atoms")
            #Extract atoms from lines
            atomsA = [(a[30:38],a[38:46],a[46:54]) for a in chains[3]]
            atomsB = [(a[30:38],a[38:46],a[46:54]) for a in chains[i]]
            atomsA = np.asarray(atomsA,dtype="float64")
            atomsB = np.asarray(atomsB,dtype="float64")
            print("\tFinding chain " + str(i) + " center")
            #Find the center of the chains
            centerA = np.sum(atomsA,axis=0)/len(atomsA)
            centerB = np.sum(atomsB,axis=0)/len(atomsB)
            print("\tCentering chain " + str(i) + "")
            A = atomsA-np.hstack(centerA)
            B = atomsB-np.hstack(centerB)
            A = np.transpose(A)
            B = np.transpose(B)
            print("\tCalculating rotation matrix for chain " + str(i) + "")
            #Calculate the H matrix
            H = np.matmul(B,np.transpose(A))
            #Decompose the H matrix by SDV
            SDV = np.linalg.svd(H)
            U = SDV[0]
            S = SDV[1]
            V = SDV[2]
            #Obtain the rotation matrix by UV^T
            subR = np.transpose(np.matmul(U,V))
            R += [subR]
            
            
            print("\tCalculating translation vector for chain " + str(i) + " with gradient descent")
            subT = np.asarray([0.0,0.0,0.0])
            error = 1000000;
            for j in range(200):
                diffCoords = [np.matmul(subR,x-centerB)+centerB+subT-y for x,y in zip(atomsB,atomsA)]
                diffCoordsX = [np.matmul(subR,x-centerB)+centerB++subT+np.asarray([0.00001,0,0])-y for x,y in zip(atomsB,atomsA)]
                diffCoordsY = [np.matmul(subR,x-centerB)+centerB++subT+np.asarray([0,0.00001,0])-y for x,y in zip(atomsB,atomsA)]
                diffCoordsZ = [np.matmul(subR,x-centerB)+centerB++subT+np.asarray([0,0,0.00001])-y for x,y in zip(atomsB,atomsA)]
                prevError = error
                error = sum(np.dot(i,i) for i in diffCoords)
                errorD = (sum(np.dot(i,i) for i in diffCoordsX)-error, sum(np.dot(i,i) for i in diffCoordsY)-error,
                sum(np.dot(i,i) for i in diffCoordsZ)-error)
                errorD = np.asarray(errorD)/0.0001
                step = 0.001 if j < 100 else 0.0000001
                subT = subT - step*errorD
                if (np.abs(error - prevError) < 0.00000000000001):
                    break;
            T += [subT]
            C += [centerB-pdbOrigin]
            print("\tFinished chain " + str(i) + "")

    #This algorithm detects what chain does the conformation belongs to and selects the aproppriate rotation matrix to apply
    #This algorithm is dependant on the Shaker.pdb chain structure. 
    #Each chain takes one quadrant of the xy-plane. This is, in counter-clockwise order:
    #Chain C, Chain B Chain D, Chain A
    chainsPositions = [2,1,3,0]

    confFile = open(confPath, "r")
    conformationLines = confFile.read()
    confFile.close();
    conformations = [[y for y in x.split("\n") if y and "#E:" not in y and "#I:" not in y] for x in conformationLines.split("\n================================================================================\n") if x]
    pdbqtLigandName = confPath.split("/")[-1].split(".")[0]

    if True:
        #Get all conformations center
        allAtoms = []
        for i in range(len(conformations)):
            conformation = conformations[i]
            allAtoms += [np.asarray((a[30:38],a[38:46],a[46:54]),dtype="float64") for a in conformation if "ATOM" in a]
        confOrigin = np.sum(allAtoms,axis=0)/len(allAtoms)
        
        #Transform all configs to align with the D chain
        for i in range(len(conformations)):
            conformation = conformations[i]
            #Get all atoms in current conformation into array
            confAtom = [np.asarray((a[30:38],a[38:46],a[46:54]),dtype="float64") for a in conformation if "ATOM" in a]
            #Find the chain this conformation belongs to by looking at the angle it makes with the x-axis
            confCenter = np.sum(confAtom,axis=0)/len(confAtom)
            #confCenter -= confOrigin
            angle = np.arctan2(confCenter[1],confCenter[0]) * 180/3.14159235657989323846
            angle = angle if angle > 0 else angle+360
            chainInfo = chainsPositions[int(angle / 90)]
            if (chainInfo != 3):
                currentR = R[chainInfo]
                currentC = C[chainInfo]
                currentT = T[chainInfo]
                confAtom = [(np.matmul(currentR, x-currentC)+currentC+currentT) for x in confAtom]
            for j in range(len(conformation)):
                if "ATOM" in conformation[j]:
                    conformation[j] = conformation[j][0:30] + "{:>8.3f}".format(confAtom[j][0]) + "{:>8.3f}".format(confAtom[j][1]) + "{:>8.3f}".format(confAtom[j][2]) + conformation[j][54:]
                conformations[i] = conformation
    confWriter = open(pdbqtLigandName+"_T.conf", "w")

    for i in range(len(conformations)):
        confWriter.write("\n".join(conformations[i]))
        confWriter.write("\n")
        confWriter.write("================================================================================")
        confWriter.write("\n")

if __name__ == '__main__':
    MashConfigurations(sys.argv[1:])
    input("Press enter to exit")