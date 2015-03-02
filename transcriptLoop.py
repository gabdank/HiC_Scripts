# parse transcripts ensemble file and then use the TSS and TES to define the interesting
# regions to count interactions from
# record the numbers of interactions and create a distribution of numbers
# repeat the same for randomly selected sites in the genome - with the same lengths of the original transcripts
# may be store as well the relative location of the random transcripts

chrLengthDict = {"chrI":15072423, "chrII":15279345,"chrIII":13783700,"chrIV":17493793,"chrV":20924149,"chrX":17718866,"chrM":13794}
chrCumulativeLengthsDictionary = {"chrI":0, "chrII":15072423,"chrIII":30351768,"chrIV":44135468,"chrV":61629261,"chrX":82553410,"chrM":100272276}

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

dict = {}
detectionsFileReadIn("/home/gabdank/Documents/January28/GLP_AVA/deduped.filtered.detections",dict)

#print len(dict)

'''for x in dict:
    print x
    print dict[x]
    break
'''
for transcript in transcriptDictionary:
    print transcript
    print transcriptDictionary[transcript]
    break