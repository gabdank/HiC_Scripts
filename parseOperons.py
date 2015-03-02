import random

chrLengthDict = {"chrI":15072423, "chrII":15279345,"chrIII":13783700,"chrIV":17493793,"chrV":20924149,"chrX":17718866,"chrM":13794, "chrMtDNA":13794}
chrCumulativeLengthsDictionary = {"chrI":0, "chrII":15072423,"chrIII":30351768,"chrIV":44135468,"chrV":61629261,"chrX":82553410,"chrM":100272276, "chrMtDNA":100272276}

def detectionsFileReadIn(fileName, dictionary):
    f = open(fileName,"r")
    mone = 0
    for l in f:
        mone = mone+1
        if (mone%100000==0):
            print "Processed "+str(mone)+" of lines in detections file " + str(fileName)
        arr = l.strip().split()
        chrA = arr[-4]
        chrB = arr[-2]

        if chrA==chrB:
            posA = int(arr[-3])
            posB = int(arr[-1])

            if posA<posB:
                left = posA
                right = posB
            else:
                left = posB
                right = posA

            key = (chrA,left)

            if not key in dictionary:
                dictionary[key]= (right, 1)
            else:
                (a,b) = dictionary[key]
                dictionary[key] = (a,b+1)

    f.close()


transcriptsDict = {}
detF = open("/home/gabdank/Documents/January28/transcripts/ensembel.transcripts","r")
for l in detF:    
    arr = l.strip().split()
    transcriptName = arr[0]
    chromosoma = "chr"+arr[-3]
    start = int(arr[-2])
    end = int(arr[-1])
    if not transcriptName in transcriptsDict:
        transcriptsDict[transcriptName]=(chromosoma, start,end)
        
detF.close()
print "The length of transcripts dictionary is:"+str(len(transcriptsDict))


histo = {}
operonsDict = {}
opF = open("operons.list","r")
outOper = open("operonCoords.list","w")
for l in opF:    
    #print l
    operonStart = 100000000
    operonEnd = -1
    
    arr = l.strip().split()
        
    flag = False
    
    for x in arr:
        if x in transcriptsDict:
            flag = True
            if operonStart>transcriptsDict[x][1]:
                operonStart=transcriptsDict[x][1]
            if operonEnd<transcriptsDict[x][2]:
                operonEnd=transcriptsDict[x][2]
            operonChr = transcriptsDict[x][0]
    if flag==True:
        key = (operonChr, operonStart,operonEnd)            
        operonsDict[key]=1    
        operonLength = operonEnd-operonStart+1
        binned = operonLength/100
        if not binned in histo:
            histo[binned]=0
        histo[binned]+=1
        if operonLength>4000:        
            outOper.write(operonChr+"\t"+str(operonStart)+"\t"+str(operonEnd)+"\n")

opF.close()
outOper.close()

print "number of operons is "+str(len(operonsDict))


randomPseudoTranscripts = {}
for name in operonsDict.keys():
    (chr,start,end) = name
    delta = end-start+1
    upperL = chrLengthDict[chr]-delta
    randoStart = random.randint(1,upperL)
    randoEnd = randoStart+delta -1
    randomPseudoTranscripts[name]=(chr,randoStart, randoEnd)




dict = {}
detectionsFileReadIn("/home/gabdank/Documents/January28/N2_DPN/deduped.filtered.detections",dict)


flankRange = 1000
histogram = {}

#for x in operonsDict.keys():
#    (chromo,start,end)=x
for x in randomPseudoTranscripts:
    (chromo,start,end) = randomPseudoTranscripts[x]
    leftStart = start-flankRange
    leftEnd = start+flankRange+1
    rightStart = end-flankRange
    rightEnd = end+flankRange+1

    mone = 0
    for x in range(leftStart,leftEnd):
        k = (chromo,x)
        if k in dict:
            score = dict[k][1]
            mone += score

    if not mone in histogram:
        histogram[mone]=1
    else:
        histogram[mone]+=1

outputF = open("/home/gabdank/Documents/January28/histogra.operon.rand","w")
for juncNum in sorted(histogram.keys()):
    outputF.write(str(juncNum)+"\t"+str(histogram[juncNum])+"\n")
outputF.close()