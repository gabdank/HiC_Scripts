import random as rnd
import numpy as np

chr2Num = {'I':1,'II':2,'III':3,'IV':4,'V':5,'X':6,'M':7,'MtDNA':7,'chrI':1, 'chrII':2, 'chrIII':3, 'chrIV':4, 'chrV':5, 'chrX':6, 'chrMtDNA':7, 'chrM':7}
num2Chr = {1:'chrI', 2:'chrII', 3:'chrIII', 4:'chrIV', 5:'chrV', 6:'chrX', 7:'chrM'}

# step 1 create 13 arrays: for each chromosome, for all except one chromosome, for all chromosomes (total of 13 arrays)
chrI_dict = []
chrII_dict = []
chrIII_dict = []
chrIV_dict = []
chrV_dict = []
chrX_dict = []
chrM_dict = []

except_chrI_dict = []
except_chrII_dict = []
except_chrIII_dict = []
except_chrIV_dict = []
except_chrV_dict = []
except_chrX_dict = []
except_chrM_dict = []
allSites_dict = []

restrDict = {-7:except_chrM_dict,-1:except_chrI_dict,-2:except_chrII_dict,-3:except_chrIII_dict,-4:except_chrIV_dict,-5:except_chrV_dict,-6:except_chrX_dict,1:chrI_dict,
             2:chrII_dict,3:chrIII_dict,4:chrIV_dict,5:chrV_dict,6:chrX_dict,7:chrM_dict,0:allSites_dict}



dpnSitesFile = open("/media/gabdank/Backup/NextSeq/MonteCarlo/avaII.sorted.bed","r")

for li in dpnSitesFile:
    arr = li.strip().split()
    chromo = arr[0]
    start = int(arr[1])
    end = int(arr[2])
    length = end-start
    position = start+(length/2)
    chrIndex = chr2Num[chromo]
    pair = (chrIndex, position)
    restrDict[0].append(pair)
    restrDict[chrIndex].append(pair)
    for k in restrDict.keys():
        if (k<0 and k!=-(chrIndex)):
             restrDict[k].append(pair)

dpnSitesFile.close()

# generate junctions list of a given length
intraChromasomalProbability = 0.999999
interchromasomalProbability = 1.0-intraChromasomalProbability


junctionsList = []

for x in range(0,10000000):
    if rnd.random()<intraChromasomalProbability: #intra chromasomal junction
        firstEnd = rnd.choice(restrDict[1])
        firstEndChr = firstEnd[0]
        secondEnd = rnd.choice(restrDict[firstEndChr])
    #else: #interchromasomal junction
    #    firstEnd = rnd.choice(restrDict[0])
    #    firstEndChr = firstEnd[0]
    #    secondEnd = rnd.choice(restrDict[-firstEndChr])
    junction = (firstEnd, secondEnd)
    junctionsList.append(junction)
    if x%1000000==0:
        print "created "+str(x) + " random junctions"
print "---------------------------"
#print len(junctionsList)

# generate matrix of a given chromosome with a given resolution
# creating matrix of the genome
chrLengthDict = {"chrI":15072423, "chrII":15279345,"chrIII":13783700,"chrIV":17493793,"chrV":20924149,"chrX":17718866,"chrM":13794}
chrCumulativeLengthsDictionary = {"chrI":0, "chrII":15072423,"chrIII":30351768,"chrIV":44135468,"chrV":61629261,"chrX":82553410,"chrM":100272276}

resolution = 50000
matrixKeyDict = {}
mone = 0
#for x in range(1,8):
x = 1
currentChromosomeLength = chrLengthDict[num2Chr[x]]
for y in range (1, currentChromosomeLength, resolution):
    #print y
    matrixKeyDict[(x,y)]=mone
    mone += 1

print len((matrixKeyDict.keys()))

# initiating 2D array
junctionsMatrix = np.zeros((len((matrixKeyDict.keys())),len((matrixKeyDict.keys()))))
#print junctionsMatrix

mone = 0
for junction in junctionsList:
    mone += 1
    if mone % 1000000 == 0 :
        print "Added "+str(mone)+" junctions to the matrix"
    #print junction
    left = junction[0]
    right = junction[1]

    leftAddress = (((left[1]/resolution)*resolution)+1)
    leftMatrixKey = (left[0],leftAddress)

    rightAddress = (((right[1]/resolution)*resolution)+1)
    rightMatrixKey = (right[0],rightAddress)

    junctionsMatrix[matrixKeyDict[leftMatrixKey], matrixKeyDict[rightMatrixKey]] += 1
    if (leftMatrixKey!=rightMatrixKey):
        junctionsMatrix[matrixKeyDict[rightMatrixKey], matrixKeyDict[leftMatrixKey]] += 1
    #break
print "Finshed to add all the junctions"
np.savetxt('/media/gabdank/Backup/NextSeq/MonteCarlo/mat.mat',junctionsMatrix,delimiter="\t")