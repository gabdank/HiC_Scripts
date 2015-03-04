# parse transcripts ensemble file and then use the TSS and TES to define the interesting
# regions to count interactions from
# record the numbers of interactions and create a distribution of numbers
# repeat the same for randomly selected sites in the genome - with the same lengths of the original transcripts
# may be store as well the relative location of the random transcripts
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




transcriptFile = open("/home/gabdank/Documents/January28/transcripts/ensembel.transcripts","r")
transcriptDictionary = {}
for l in transcriptFile:
    transcriptArray = l.strip().split()
    transcriptName = transcriptArray[0]
    transcriptChromosome = "chr"+transcriptArray[2]
    transcriptStart = int(transcriptArray[3])
    transcriptEnd = int(transcriptArray[4])
    if not transcriptName in transcriptDictionary:
        transcriptDictionary[transcriptName] = (transcriptChromosome,transcriptStart,transcriptEnd)
transcriptFile.close()
print "Number of transcripts was:"+str(len(transcriptDictionary))

#Looking into length distribution of transcripts
transcriptLengths = {}
binSize = 25
for name in transcriptDictionary:
    (chr,start,end) = transcriptDictionary[name]
    delta = end-start+1
    binnedDelta = delta/binSize
    if not binnedDelta in transcriptLengths:
        transcriptLengths[binnedDelta]=1
    else:
        transcriptLengths[binnedDelta]+=1


randomPseudoTranscripts = {}
for name in transcriptDictionary:
    (chr,start,end) = transcriptDictionary[name]
    delta = end-start+1
    upperL = chrLengthDict[chr]-delta
    randoStart = random.randint(1,upperL)
    randoEnd = randoStart+delta -1
    randomPseudoTranscripts[name]=(chr,randoStart, randoEnd)


# getting to the detections counts

dict = {}
detectionsFileReadIn("/home/gabdank/Documents/January28/GLP_AVA/deduped.filtered.detections",dict)
detectionsFileReadIn("/home/gabdank/Documents/January28/GLP_DPN/deduped.filtered.detections",dict)

#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_FED_AVA_1/deduped.filtered.detections",dict)
#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_FED_AVA_2/deduped.filtered.detections",dict)
#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_FED_DPN_1/deduped.filtered.detections",dict)
#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_FED_DPN_2/deduped.filtered.detections",dict)

#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_STARVED_AVA_1/deduped.filtered.detections",dict)
#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_STARVED_AVA_2/deduped.filtered.detections",dict)
#detectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_STARVED_DPN_1/deduped.filtered.detections",dict)
#etectionsFileReadIn("/media/gabdank/Disk3/HiC/pipeline/HISEQ_STARVED_DPN_2/deduped.filtered.detections",dict)


flankRange = 1000
histogram = {}


for transcript in transcriptDictionary:
    (chromo,start,end) = transcriptDictionary[transcript]
#for transcript in randomPseudoTranscripts:
#    (chromo,start,end) = randomPseudoTranscripts[transcript]

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

outputF = open("/home/gabdank/Documents/January28/histogra.glp.ava.dpn","w")
for juncNum in sorted(histogram.keys()):
    outputF.write(str(juncNum)+"\t"+str(histogram[juncNum])+"\n")
outputF.close()

