def show_exception_and_exit(exc_type, exc_value, tb):
    import traceback
    traceback.print_exception(exc_type, exc_value, tb)
    input("Press key to exit.")
    sys.exit(-1)

import sys
import numpy as np
import MashConfigs
import os
import os.path as path

from scipy.stats import norm
import operator
import matplotlib.pyplot as plt
import itertools

sys.excepthook = show_exception_and_exit

if (len(sys.argv) != 3):
    raise ValueError("Specified more or less than 2 files. Required 2 files. %.conf and %.pdb files")
    
#Set docking and pdb file paths
confPath = [x.replace("\\","/") for x in sys.argv if ".conf" in x][0]
pdbPath = [x.replace("\\","/") for x in sys.argv if ".pdb" in x][0]

if (confPath == "" or pdbPath == ""):
    raise ValueError("Missing one file")

#Mash all configurations together
MashConfigs.MashConfigurations([confPath, pdbPath])


confFile = open(confPath, "r")
conformationLines = confFile.read()
dockData = [[y.split(":")[-1] for y in x.split("\n") if "#" in y] for x in conformationLines.split("\n================================================================================\n") if x]

print("Reading docking energy and inhibition data")
#Create data structure to hold #conformation, energy and inhibition values
rawData = [(idx,float(x[0]),float(x[1])) for idx,x in enumerate(dockData)]
fixedData = [(a[0],0.592089376625*np.log(a[2]/1000),a[2]) for a in rawData]
sortedData = sorted(fixedData, key=operator.itemgetter(1),reverse=False)

print("Reading transformed conformations")
#Read transformed conformations
transformedConfFile = confPath.replace(".conf", "_T.conf")
tconfFile = open(transformedConfFile, "r")
conformationLines = tconfFile.read()
tconfFile.close();
conformations = [[y for y in x.split("\n") if "#E:" not in y and "#I:" not in y and y] \
for x in conformationLines.split("\n================================================================================\n") if x]
#Get minimun energy conformation for RMSE reference
minEnergyConfNumber = sortedData[0][0]
minEnergyConf = conformations[minEnergyConfNumber]
minConfAtom = [np.asarray((a[30:38],a[38:46],a[46:54]),dtype="float64") for a in minEnergyConf if "ATOM" in a]

print("Calculating conformations RMSE")
#Calculate RMSE clusters
rmse = []
for confIndex in range(len(conformations)):
    conformation = conformations[confIndex]
    #Get all atoms in current conformation into array
    confAtom = [np.asarray((a[30:38],a[38:46],a[46:54]),dtype="float64") for a in conformation if "ATOM" in a]
    difference = [np.dot((x-y),(x-y)) for x,y in zip(minConfAtom, confAtom)]
    partialRMSE = np.sqrt(sum(difference)/len(difference))
    rmse += [(confIndex,partialRMSE)]

sortedRMSE = sorted(rmse, key=operator.itemgetter(1),reverse=False)
rmseGaps = []
for i in range(len(rmse)-1):
    rmseGaps += [(i,sortedRMSE[i+1][1] - sortedRMSE[i][1])]
sortedRMSEGaps = sorted(rmseGaps, key=operator.itemgetter(1),reverse = True)

clusterTest = []
for i in range(len(sortedRMSEGaps)):
    clusterValue = sortedRMSE[sortedRMSEGaps[i][0]+1][1]/sortedRMSEGaps[i][1]
    clusterTest += [clusterValue]
lastGap = -1
for i in range(len(clusterTest)):
    if (clusterTest[i] <= 2):
        lastGap = i
    else:
        break

rmseConfClusters = []   
if (lastGap != -1):
    rmseBreaks = sorted(sortedRMSEGaps[:lastGap+1],key=operator.itemgetter(0))
    rmseBreaksStart = 0;
    for i in range(len(rmseBreaks)):
        rmseConfClusters += [sortedRMSE[rmseBreaksStart:rmseBreaks[i][0]+1]]
        rmseBreaksStart = rmseBreaks[i][0]+1
    rmseConfClusters += [sortedRMSE[rmseBreaksStart:]]
else:
    rmseConfClusters += [sortedRMSE]


dataClusters = []
for i in range(len(rmseConfClusters)):
    temp = []
    for j in range(len(rmseConfClusters[i])):
        confID = rmseConfClusters[i][j][0]
        temp += [(fixedData[confID][1],fixedData[confID][2],rmse[confID][1])]
    dataClusters += [(str(i),temp)]

#Obtain all combinations of the clusters and add them to the cluster list
clusterIndexList = list(range(len(dataClusters)))
LIndexList = list(range(2,len(dataClusters)+1))
for L in LIndexList:
    for subset in itertools.combinations(clusterIndexList, L):
        dataTemp = []
        indexTemp = ""
        for i in range(len(subset)):
            dataTemp += dataClusters[subset[i]][1]
            indexTemp += dataClusters[subset[i]][0]+"+"
        indexTemp = indexTemp[:-1]
        dataClusters += [(indexTemp, dataTemp)]


clustersInfo = ""
for i in range(len(dataClusters)):
    clustersInfo += "Cluster "+dataClusters[i][0]+" data:" + "\n"
    print("Cluster "+dataClusters[i][0]+" data:")
    mu = np.mean([x[0] for x in dataClusters[i][1]])
    std = np.std([x[0] for x in dataClusters[i][1]])
    clustersInfo += "\tEnergy mean: " + str(mu) + "\n"
    clustersInfo += "\tEnergy error: " + str(std) + "\n"
    print ("\tEnergy mean: " + str(mu))
    print ("\tEnergy error: " + str(std))
    
    mu = np.mean([x[1] for x in dataClusters[i][1]])
    std = np.std([x[1] for x in dataClusters[i][1]])
    clustersInfo += "\tInhibition mean: " + str(mu) + "\n"
    clustersInfo += "\tInhibition error: " + str(std) + "\n"
    print ("\tInhibition mean: " + str(mu))
    print ("\tInhibition error: " + str(std))
    
    mu = np.mean([x[2] for x in dataClusters[i][1]])
    std = np.std([x[2] for x in dataClusters[i][1]])
    clustersInfo += "\tRMSE mean: " + str(mu) + "\n"
    clustersInfo += "\tRMSE error: " + str(std) + "\n"
    print ("\tRMSE mean: " + str(mu))
    print ("\tRMSE error: " + str(std))
    
    clustersInfo += "\n"
    print("\n")

clustersInfoFile = open("./ClusterInfo.txt","w")
clustersInfoFile.write(clustersInfo)
clustersInfoFile.close()

pdbqtLigandName = confPath.split("/")[-1].split(".")[0]

confClusterFolders = "./"+pdbqtLigandName+"_ConfFiles"
if (not os.path.exists(confClusterFolders)):
    os.mkdir(confClusterFolders)
for i in range(len(rmseConfClusters)):
    currentClusterFolder = confClusterFolders+"/Cluster_"+str(i)
    os.mkdir(currentClusterFolder)
    for j in range(len(rmseConfClusters[i])):
        confWriter = open(currentClusterFolder + "/Cluster" + str(rmseConfClusters[i][j][0]) + ".pdb" , "w")
        confWriter.write("\n".join(conformations[rmseConfClusters[i][j][0]]))
        confWriter.write("\n")
        confWriter.write("================================================================================")
        confWriter.write("\n")
    


input("Press enter to exit")