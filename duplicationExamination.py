import sys

'''
Function processDetectionsFile takes as parameter the path to the file containing the detections (currently not filtered
Detections would contain chunks from the reference genome aling with junctions.
Each detection consists of chrA posA strandA and same for the other side of the detection (chrB, posB, strandB)
We would convert each side into ABSOLUTE numerical value (using the chromosome lengths etc. and then we will make a pair in which we will
store at the first place the smaller ABSOLUTE value, and the sign of th evalue will be dictated by the strandidness
Pair of two numbers will be the key to the dictionary of different detections.
The dictionary will store the number of occurrences of such a "pair".

'''

def processDetectionsFile(fileName):
    dict2Return = {}
    handle = open(fileName,"r")
    for lin in handle:
        print lin
        break
    handle.close()
    return dict2Return


if len(sys.argv)<2:
    print ""
    print "ERROR in the input to the python script call"
    print "ERROR, please use the path to the detections file as first parameter to this function"
    print "ERROR in the input to the python script call"
    print ""
else:
    detectionsFileName = sys.argv[1]
    processDetectionsFile(detectionsFileName)
