import numpy as np
oneMatrixFile = open("/media/gabdank/Backup/NextSeq/AF_SOL_597/N2_DPN/n2.dpn.50000.chrI.norm.noUpperLine.noFirstTwoCols","r")
twoMatrixFile = open("/media/gabdank/Backup/NextSeq/AF_SOL_597/GLP_DPN/glp.dpn.50000.chrI.norm.noUpperLine.noFirstTwoCols","r")
hugeL = []
for l in oneMatrixFile:
    arr = l.strip().split()
    hugeL.append(arr)
one = np.array(hugeL)

hugeL = []
for l in twoMatrixFile:
    arr = l.strip().split()
    hugeL.append(arr)
two = np.array(hugeL)

oneMatrixFile.close()
twoMatrixFile.close()


mixMat = open("/media/gabdank/Backup/NextSeq/AF_SOL_597/n2.glp.50000.chrI","w")
for x in range(len(one)):
    for y in range(len(one)):
        if (x<y):
            mixMat.write(str(one[x,y])+"\t")
        else:
            if (x==y):
                mixMat.write("0\t")
            else:
                mixMat.write(str(two[x,y])+"\t")

    mixMat.write("\n")
mixMat.close()

