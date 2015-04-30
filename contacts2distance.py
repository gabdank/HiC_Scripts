'''
Script allows to extract distribution of junction numbers, depending on the distance
between "bins" at the genome.
Normalization is to the totla number of junctions in a given experiment.
'''
import numpy as np

'''
# Starting from matrix of ONE chromosome
'''

def median(lst):
    return np.median(np.array(lst))

def singleChromosomeMatrixToChart(binSize, valuesMatrix, chartFileName):
    distanceDict = {}
    rowIndex = 0
    colIndex = 0
    length = len(valuesMatrix[0])
    totalNumberOfJunctions = 0
    while (rowIndex<length):
        colIndex = rowIndex
        distance = binSize
        while (colIndex<length):
            if not distance in distanceDict:
                distanceDict[distance]=[]        
            distanceDict[distance].append(valuesMatrix[rowIndex,colIndex])
            totalNumberOfJunctions += valuesMatrix[rowIndex,colIndex]
            colIndex +=1
            distance += binSize        
        rowIndex +=1

    f = open(chartFileName,"w")
    for k in sorted(distanceDict.keys()):
        f.write(str(k) +"\t" + str(median(distanceDict[k])/float(totalNumberOfJunctions))+"\n")
    f.close()
    

def getCoordinatesFromLine(line):
    arr = line.strip().split()
    return (arr[-4],arr[-3],arr[-2],arr[-1])
    
def sameChromosome(coordinatesTuple):
    if coordinatesTuple[0]==coordinatesTuple[2]:
        return True
    return False
    
def getDistance(coordinatesTuple):
    return abs(int(coordinatesTuple[1])-int(coordinatesTuple[3]))
    

def junctionsFileToChart(filteredJunctionsFileName,chartFileName,binSize):
    junctionsFile = open(filteredJunctionsFileName,"r")
    outputFile = open(chartFileName,"w")
    junctionsCounter = 0    
    dictionaryOfDistances = {}
    for l in junctionsFile:
        junctionsCounter += 1
        coordinatesTuple = getCoordinatesFromLine(l)
        if sameChromosome(coordinatesTuple)==True:        
            distance = getDistance(coordinatesTuple)    
            binnedDistance = distance/binSize
            if not binnedDistance in dictionaryOfDistances:
                dictionaryOfDistances[binnedDistance] = 1
            else:
                dictionaryOfDistances[binnedDistance] += 1
        if junctionsCounter%100000==0:
            print "Processed "+str(junctionsCounter)+" junctions"
            
    junctionsFile.close()
    
    for binnedDistance in sorted(dictionaryOfDistances.keys()):
        outputFile.write(str(binnedDistance*binSize) + "\t"+str(float(dictionaryOfDistances[binnedDistance])/float(junctionsCounter))+"\n")
    outputFile.close()
    
#junctionsFile = "GLP_AVA/deduped.filtered.detections"
#junctionsFileToChart(junctionsFile,"zopa_glp_ava",5000)

matrixFile = "/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/glp.dpn.chrIII.raw.2K.noUpperLine.noFirstTwoCols"
chartFile = "/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/glp.dpn.chrIII.raw.5K.chart"
oneMatrixFile = open(matrixFile,"r")

hugeL = []
for l in oneMatrixFile:
    arr = [float(i) for i in l.strip().split()]
    hugeL.append(arr)
one = np.array(hugeL)

oneMatrixFile.close()

singleChromosomeMatrixToChart(5000, one, chartFile)