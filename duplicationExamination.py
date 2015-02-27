import sys


chrLengthDict = {"chrI":15072423, "chrII":15279345,"chrIII":13783700,"chrIV":17493793,"chrV":20924149,"chrX":17718866,"chrM":13794}
chrCumulativeLengthsDictionary = {"chrI":0, "chrII":15072423,"chrIII":30351768,"chrIV":44135468,"chrV":61629261,"chrX":82553410,"chrM":100272276}


'''
Function processDetectionsFile takes as parameter the path to the file containing the detections (currently not filtered
Detections would contain chunks from the reference genome aling with junctions.
Each detection consists of chrA posA strandA and same for the other side of the detection (chrB, posB, strandB)
We would convert each side into ABSOLUTE numerical value (using the chromosome lengths etc. and then we will make a pair in which we will
store at the first place the smaller ABSOLUTE value, and the sign of th evalue will be dictated by the strandidness
Pair of two numbers will be the key to the dictionary of different detections.
The dictionary will store the number of occurrences of such a "pair".
M03109:43:000000000-ACGWU:1:1101:10000:10610	chrV	3158715	False	TAGATGAATTCCATCTGGAATCGGGTCGACC	chrV	3158588	True	CGTTTCGACAATCTCCTTCTTCATGC

'''

def processDetectionsFile(fileName):
    dict2Return = {}
    handle = open(fileName,"r")
    counter = 0
    for lin in handle:
        counter += 1
        if (counter % 1000000==0):
            print "Processed "+str(counter)+" lines in  detections file"
        arr = lin.strip().split()

        chrA = arr[1]
        posA = int(arr[2])
        strandA = arr[3]
        if strandA=='True':
            signA = 1
        else:
            signA = -1

        chrB = arr[5]
        posB = int(arr[6])
        strandB = arr[7]
        if strandB=='True':
            signB = 1
        else:
            signB = -1

        absoluteA = chrCumulativeLengthsDictionary[chrA]+posA
        absoluteB = chrCumulativeLengthsDictionary[chrB]+posB

        if absoluteA<absoluteB:
            key = (signA*absoluteA,signB*absoluteB)
        else:
            key = (signB*absoluteB,signA*absoluteA)
        if not key in dict2Return:
            dict2Return[key]=1
        else:
            dict2Return[key]+=1

    #print "Scanned "+str(counter)+" entries"
    handle.close()
    toReturn = (dict2Return, counter)
    return toReturn


if len(sys.argv)<2:
    print ""
    print "ERROR in the input to the python script call"
    print "ERROR, please use the path to the detections file as first parameter to this function"
    print "ERROR in the input to the python script call"
    print ""
else:
    detectionsFileName = sys.argv[1]
    (dictionary, totalNumber) = processDetectionsFile(detectionsFileName)
    print "Number of unique entries: "+str(len(dictionary.keys()))
    print "Nuber of total entries: "+str(totalNumber)
    sum = float(0)
    mone = 0
    for k in dictionary:
        mone +=1
        if (mone % 1000000==0):
            print "Another million used"
            break
        observed = dictionary[k]
        p_of_k = float(observed)/float(totalNumber)
        p_of_k_square = p_of_k*p_of_k
        sum += p_of_k_square
    print "The sigma of squared probabilities = "+str(sum)