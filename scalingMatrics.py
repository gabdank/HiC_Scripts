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
