# parse transcripts ensemble file and then use the TSS and TES to define the interesting
# regions to count interactions from
# record the numbers of interactions and create a distribution of numbers
# repeat the same for randomly selected sites in the genome - with the same lengths of the original transcripts
# may be store as well the relative location of the random transcripts

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

for transcript in transcriptDictionary:
    print transcript
    print transcriptDictionary[transcript]
    break