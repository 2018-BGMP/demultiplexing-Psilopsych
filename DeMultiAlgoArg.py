
# coding: utf-8

# In[29]:


#!/usr/bin/env python3

from pathlib import Path
import argparse
import gzip

def get_arguments():
        parser = argparse.ArgumentParser()
        parser.add_argument("-r1", "--read1", help="Forward read file", type = str, required = True)
        parser.add_argument("-i1", "--index1", help="Forward index file", type = str, required = True)
        parser.add_argument("-r2", "--read2", help="Reverse read file", type = str, required = True)
        parser.add_argument("-i2", "--index2", help="Reverse index file", type = str, required = True)
        parser.add_argument("-t", "--threshold", help="Specify the k-mer size you intend to work with", type = float, required = False, default = 25,        choices = range(0, 40))
        parser.add_argument("-b", "--barcodes", help="Barcode mapping file", type = str, required = False, default = "indexes.txt")
        return parser.parse_args()

def reverseComplement(dna):
    """This function returns the reverse complement of a DNA sequence with ATGC's present"""
    dna = dna.upper()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def convert_phred(letter):
    """Converts a single character in binary into a phred score. Ord function does this """
    phred = ord(letter) - 33 #Phred 33 encoding
    return phred

def indexDict(indexFile):
    """ This function returns a dictionary of iD barcode and barcode sequence
    in a dictionary with the sequence as the key, and ID as value. """
    index = {}
    with open(indexFile) as table:
        for line in table:
            if line.startswith('sample') == False:
                line = line.strip().split()
                index[line[4]] = line[3]
    return index

def writeOutFq(fRead, rRead, indexSeq, fileDict):
    """ This function is used to write out to two files, default files are those specified in parameter. The name of output file are
    found in index dictionary. N in dictionary represents bad files."""
    #Find appropiate ID from forward sequence ID
    ID = indexLib[indexSeq]
    forward = ID + "f"
    reverse = ID + "r"
    forward = fileDict[forward]
    reverse = fileDict[reverse]
    forward.write(fRead)
    reverse.write(rRead)
     

# barcodes = "indexes.txt"
barcodes = args.barcodes
#This dictionary holds barcode sequence as key and respective ID as value
indexLib = indexDict(barcodes)   
#Add N to dictionary for all bad reads 
indexLib["N"] = "unMatched"
# args = get_arguments()
# bioRead1 = args.read1
# indexRead1 = args.index1
# bioRead2 = args.read2
# indexRead2 = args.index2
# qsThreshold = args.threshold

#Test Files
bioRead1 = 'test.fastq.gz'
indexRead1 = 'testIndex.fastq.gz'
bioRead2 = 'testR2.fastq.gz'
indexRead2 = 'testIndexR2.fastq.gz'
qsThreshold = 25


#This dictionary holds ID + forward/reverse read as key, and respective opened file as values.
files = {}
for k,v in indexLib.items():
    fileF = open(v + "Forward.fastq", "a+")
    fileR = open(v + "Reverse.fastq", "a+")
    files[v + "f"] = fileF
    files[v + "r"] = fileR
 
    
with gzip.open(bioRead1, "rt") as bio, gzip.open(indexRead1, "rt") as index, gzip.open(indexRead2, "rt") as indexr, gzip.open(bioRead2, "rt") as bior:
    #The following variables keep track of matching statistics 
    recordsRead = 0
    #Matched barcodes that surpass QS threshold
    matchedInd = 0
    
    #Mean quality score of barcodes dont surpass or equal threshold value, written to bad files
    belowThres = 0
    
    #Mismatched barcodes, written to bad files
    unmatchedInd = 0
    
    # N in quality scores, automatically written to bad files
    badRead = 0
    
    #Matched indexes not found in index library
    notFound = 0
    
    
    while True:
        #Bio reads 1
        b1 = bio.readline().strip()
        b2 = bio.readline().strip()
        b3 = bio.readline().strip()
        b4 = bio.readline().strip()
        #Index reads 1
        i1 = index.readline().strip()
        i2 = index.readline().strip()
        i3 = index.readline().strip()
        i4 = index.readline().strip()
        #Index reads 2
        i5 = indexr.readline().strip()
        i6 = indexr.readline().strip()
        i7 = indexr.readline().strip()
        i8 = indexr.readline().strip()
        #Bio reads 2
        b5 = bior.readline().strip()
        b6 = bior.readline().strip()
        b7 = bior.readline().strip()
        b8 = bior.readline().strip()
        #Break loop if you reach the end of either of the files. Should(theoretically) process same amount of lines for each file
        if b1 == "" or i1 == "" or i5 == "" or b5 == "":
            break
        read = (b1 + ":" + i2) + '\n' + b2 + '\n' + b3 + '\n' + b4 + '\n'
        readR = (b5 + ":" + i6) + '\n' + b6 + '\n' + b7 + '\n' + b8 + '\n'
        recordsRead += 1
        if 'N' in i2 or 'N' in i6:
            writeOutFq(read, readR, "N", files)
            badRead += 1
            continue
        #Running sum of quality scores to calc Avg
        totalPhred = 0
        totalPhredR = 0
        sumF = 0
        sumR = 0
        #For loops used to calculate mean quality score. Only done for one since they match
        for char in i4:
            totalPhred += convert_phred(char)
            sumF += 1
        mean = totalPhred / sumF
        for char in i8:
            totalPhredR += convert_phred(char)
            sumR += 1
        meanR = totalPhredR / sumR
        if mean >= qsThreshold and meanR >= qsThreshold:
            if i2 == reverseComplement(i6):    
                if i2 in indexLib:
                    matchedInd += 1 
                    readR = (b5 + ":" + i2) + '\n' + b6 + '\n' + b7 + '\n' + b8 + '\n'
                    writeOutFq(read, readR, i2, files)
                    continue
                else:
                    notFound += 1
                    writeOutFq(read, readR, "N", files)
                    continue

            else:
                unmatchedInd += 1
                writeOutFq(read, readR, "N", files)
                continue
        else:
            belowThres += 1
            writeOutFq(read, readR, "N", files)
            continue
#This loop closes all files stores in files dictionary
for ID, file in files.items():
    file.close()
    
print("Write statisitcs")
print("Records Read: ", recordsRead)
print("Bad reads(Ns): ", badRead)
print("Below threshold: ", belowThres)
print("unmatched indices: ", unmatchedInd)
print("Not found in indexes file: ", notFound)
print("Matched reads: ", matchedInd)

