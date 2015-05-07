import numpy as np

# read in matrix and find out the maximal value
def findMaxInMatrix(matrixFileName):
    mF = open(matrixFileName,"r")
    mone = 0
    matrixValues = []
    for li in mF:
        vals = li.strip().split()
        valsList = []
        for a in vals:
            valsList.append(float(a))

        matrixValues.append(valsList)
        mone += 1
    mF.close()

    npMatrix = np.array(matrixValues)
    print "Values range from: " +str(np.amin(npMatrix))+ " to "+str(np.amax(npMatrix))
    return npMatrix

def updateMinMax(mat, minV, maxV):
    mat[mat>maxV] = maxV
    mat[mat<minV] = minV
    return mat
'''
fileName="/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrX.norm.noUpperLine.noFirstTwoCols"
chrX = findMaxInMatrix(fileName)
fileName="/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrV.norm.noUpperLine.noFirstTwoCols"
chrV = findMaxInMatrix(fileName)
fileName="/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrIV.norm.noUpperLine.noFirstTwoCols"
chrIV = findMaxInMatrix(fileName)
fileName="/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrIII.norm.noUpperLine.noFirstTwoCols"
chrIII = findMaxInMatrix(fileName)
fileName="/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrII.norm.noUpperLine.noFirstTwoCols"
chrII = findMaxInMatrix(fileName)
fileName="/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrI.norm.noUpperLine.noFirstTwoCols"
chrI  = findMaxInMatrix(fileName)


updateMinMax(chrI,-4,2)
updateMinMax(chrII,-4,2)
updateMinMax(chrIII,-4,2)
updateMinMax(chrIV,-4,2)
updateMinMax(chrV,-4,2)
updateMinMax(chrX,-4,2)



np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/chrI.out',chrI, delimiter='\t')
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/chrII.out',chrII, delimiter='\t')
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/chrIII.out',chrIII, delimiter='\t')
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/chrIV.out',chrIV, delimiter='\t')
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/chrV.out',chrV, delimiter='\t')
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/chrX.out',chrX, delimiter='\t')


combinedX = "/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrX"
combX = findMaxInMatrix(combinedX)
updateMinMax(combX,-4,2)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/comb.chrX.out',combX,delimiter="\t")

combinedV = "/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrV"
combV = findMaxInMatrix(combinedV)
updateMinMax(combV,-4,2)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/comb.chrV.out',combV,delimiter="\t")

combinedIV = "/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrIV"
combIV = findMaxInMatrix(combinedIV)
updateMinMax(combIV,-4,2)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/comb.chrIV.out',combIV,delimiter="\t")

combinedIII = "/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrIII"
combIII = findMaxInMatrix(combinedIII)
updateMinMax(combIII,-4,2)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/comb.chrIII.out',combIII,delimiter="\t")


combinedII = "/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrII"
combII = findMaxInMatrix(combinedII)
updateMinMax(combII,-4,2)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/comb.chrII.out',combII,delimiter="\t")


combinedI = "/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrI"
combI = findMaxInMatrix(combinedI)
updateMinMax(combI,-4,2)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/comb.chrI.out',combI,delimiter="\t")



genomeFile = "/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/genome.step2.txt"
genome= findMaxInMatrix(genomeFile)
updateMinMax(genome,0,1000)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/genome.step2.mod.txt',genome,delimiter="\t")
'''

simpleI = "/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.chrI.simpleNorm.50K.noUpperLine.noFirstTwoCols"
combI = findMaxInMatrix(simpleI)
updateMinMax(combI,-4,4)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.simple.chrI.out',combI,delimiter="\t")

simpleX = "/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.chrX.simpleNorm.50K.noUpperLine.noFirstTwoCols"
combX = findMaxInMatrix(simpleX)
updateMinMax(combX,-4,4)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.simple.chrX.out',combX,delimiter="\t")

simpleI = "/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/glp.dpn.chrI.simpleNorm.50K.noUpperLine.noFirstTwoCols"
combI = findMaxInMatrix(simpleI)
updateMinMax(combI,-4,4)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/glp.simple.chrI.out',combI,delimiter="\t")

simpleX = "/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/glp.dpn.chrX.simpleNorm.50K.noUpperLine.noFirstTwoCols"
combX = findMaxInMatrix(simpleX)
updateMinMax(combX,-4,4)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/glp.simple.chrX.out',combX,delimiter="\t")


genomeFile = "/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.genome.simpleNorm.100K.noUpperLine.noFirstTwoCols"
genome= findMaxInMatrix(genomeFile)
updateMinMax(genome,-5,5)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/genome.simpleNorm',genome,delimiter="\t")

genomeFile = "/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/glp.dpn.genome.simpleNorm.100K.noUpperLine.noFirstTwoCols"
genome= findMaxInMatrix(genomeFile)
updateMinMax(genome,-5,5)
np.savetxt('/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/genome.simpleNorm',genome,delimiter="\t")