def show_exception_and_exit(exc_type, exc_value, tb):
    import traceback
    traceback.print_exception(exc_type, exc_value, tb)
    input("Press key to exit.")
    sys.exit(-1)

import sys
import os
import numpy as np
import re

from scipy.stats import norm
import operator
sys.excepthook = show_exception_and_exit

def merge(list1, list2):
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list

#if (len(sys.argv) != 3):
#    raise ValueError("Specified more or less than 2 files. Required 2 files. %.dlg and %.pdb files with CONECT records")
   
#Set docking and pdb file paths
dlgPath = [x.replace("\\","/") for x in sys.argv if ".dlg" in x]
dlgPath = "" if (len(dlgPath) == 0) else dlgPath[0]
pdbPath = [x.replace("\\","/") for x in sys.argv if ".pdb" in x]
pdbPath = "" if (len(pdbPath) == 0) else pdbPath[0]

if (dlgPath == ""):
    raise ValueError("Missing .dlg file")

if (pdbPath == ""):
    print("No CONECT records, will generate file without. Visualization on software might be compromised")
    conectRegisters = ""
else:
#Extract pdb CONECT registers
    pdbLines = []
    pdbFile = open(pdbPath,"r")
    pdbLines += pdbFile.read().split("\n")
    pdbFile.close()
    conectLines = [x for x in pdbLines if "CONECT" in x]
    conectRegisters = "\n".join(conectLines)



#Extract dockings from dock file
lines = []
with open(dlgPath,"r") as file:
    lines += file.read().split("\n")


#Extract ligand file name
pdbqtLine = [x for x in lines if re.search("move.+\.pdbqt", x)][0]
pdbqtLigandName = re.findall("(?<=move ).*(?=.pdbqt)",pdbqtLine)[0]

if (pdbqtLigandName == ""):
    raise AssertionError(".pdbqt name not found. Check .dlg format.")
print(f".pdbqt file found. File name is {pdbqtLigandName}")

#Extract all the energy records and inhibition records
energyRegex = re.compile("[\d.-]*(?= kcal\/mol)")
inhibitionRegex = re.compile("[\d.-]*(?= mM)|[\d.-]*(?= uM)")
energy = [float(re.findall(energyRegex,x)[0])*1000 for x in lines if "Estimated Free Energy of Binding" in x and "DOCKED" in x]
inhibition = [float(re.findall(inhibitionRegex,x)[0])*(0.001 if ('mM' in x) else (0.000001 if ('uM' in x) else 1)) for x in lines if "Estimated Inhibition Constant, Ki" in x and "DOCKED" in x]

if len(inhibition) == 0:
    print("There were no inhibition records found. Will only print energy records.")
elif len(energy) == len(inhibition):
    print("Energy and inhibition records count match. Proceed")
else:
    raise AssertionError("Record count of energy and inhibition do not match.")
print("\tEnergy records found: ", len(energy))
print("\tInhibition records found: ", len(inhibition))


ligandAtomsNum = len([x for x in lines if "INPUT-LIGAND-PDBQT: ATOM" in x])
mergeFactor = len([x for x in lines if "INPUT-LIGAND-PDBQT: TORSDOF" in x])
ligandAtomsNum /= mergeFactor; #If multiple dlgs were merged, you need to remove copies of the "INPUT-..." definition
ligandAtomsNum = int(ligandAtomsNum+1) #Add 1 to account for the TER atom
print("Atoms in this ligand: ", ligandAtomsNum)

conformations = []
nConf = [x.replace("DOCKED: ","") for x in lines if "DOCKED: ATOM" in x or "DOCKED: TER" in x]
dockingNum = int(len(nConf)/ligandAtomsNum)
for i in range(dockingNum):
    conformations += [nConf[i*ligandAtomsNum:i*ligandAtomsNum + ligandAtomsNum]]

confWriter = open(pdbqtLigandName+".conf", "w")
for i in range(len(conformations)):
    confWriter.write("#E:"+str(energy[i])+"\n")
    if (len(inhibition) != 0):
        confWriter.write("#I:"+str(inhibition[i])+"\n")
    confWriter.write("\n".join(conformations[i]))
    confWriter.write("\n")
    confWriter.write(conectRegisters)
    confWriter.write("\n")
    confWriter.write("================================================================================")
    confWriter.write("\n")


#Stop program from exiting
input("Press enter to exit")