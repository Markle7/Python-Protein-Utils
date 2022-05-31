def show_exception_and_exit(exc_type, exc_value, tb):
    import traceback
    traceback.print_exception(exc_type, exc_value, tb)
    input("Press key to exit.")
    sys.exit(-1)


import sys
import os
import numpy as np
import re

confPath = [x.replace("\\","/") for x in sys.argv if ".conf" in x]
confPath = "" if (len(confPath) == 0) else confPath[0]

if (confPath == ""):
    raise ValueError("Missing .conf file")

lines = []
with open(confPath, "r") as file:
    lines += file.read().split("\n")

energyRegex = re.compile('(?<=#E:).*')
inhibitionRegex = re.compile('(?<=#I:).*')

energy = [float(re.findall(energyRegex,x)[0]) for x in lines if "#E:" in x]
inhibition = np.array([float(re.findall(inhibitionRegex,x)[0]) for x in lines if "#I:" in x])
if(len(inhibition) != 0):
    energy = np.array([0.592089376625*np.log(x) for x in inhibition])
inhibition *= 1000

energyMean = np.mean(energy)
energySTD = np.std(energy)
inhibitionMean = np.mean(inhibition)
inhibitionSTD = np.std(inhibition)

print("Energy mean: ",energyMean, "kcal/mol  |  Energy std: ", energySTD," kcal/mol");
print("Inhibition mean: ",inhibitionMean, "mM  |  Inhibition std: ", inhibitionSTD," mM");

#Stop program from exiting
input("Press enter to exit")