
# coding: utf-8

# In[22]:


from pathlib import Path


def convert_phred(letter):
    """Converts a single character into a phred score. Ord function does this """
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


def writeOutFq(fRead, rRead, forwardFile = 'badForward.fastq', reverseFile = 'badReverse.fastq'):
    """ This function is used to write out to two files, default files are those specified in parameter"""
    with open(forwardFile, 'a+') as forward, open(reverseFile, 'a+') as reverse:
        forward.write(fRead)
        reverse.write(rRead)
    
indexLib = indexDict("indexes.txt")   
    
bioRead1 = 'test.fastq'
indexRead1 = 'testIndex.fastq'
bioRead2 = 'testR2.fastq'
indexRead2 = 'testIndexR2.fastq'
qsThreshold = 25
with open(bioRead1) as bio, open(indexRead1) as index, open(indexRead2) as indexr, open(bioRead2) as bior:
    #The following variables keep track of matching statistics
    
    #Matched barcodes that surpass QS threshold
    matchedInd = 0
    
    #Mean quality score of barcodes dont surpass or equal threshold value, written to bad files
    belowThres = 0
    
    #Mismatched barcodes, written to bad files
    unmatchedInd = 0
    
    # N in quality scores, automatically written to bad files
    badRead = 0
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
        if b1 == "":
            break
        read = b1 + '\n' + b2 + '\n' + b3 + '\n' + b4 + '\n'
        readR = b5 + '\n' + b6 + '\n' + b7 + '\n' + b8 + '\n'
        print(readR)
        if 'N' in i2 or 'N' in i6:
            print('N in barcode, writing to bad file')
            writeOutFq(read, readR)
            badRead += 1
            continue  
        if i2 == i6:
            #Running sum of quality scores
            totalPhred = 0
            chars = 0
            #For loops used to calculate mean quality score. Only done for one since they match
            for char in i2:
                totalPhred += convert_phred(char)
                chars += 1
            mean = totalPhred / chars
            if mean >= qsThreshold:
                matchedInd += 1
                print("Matched barcode")
                #Find barcode in dictictionary 
                id = indexLib[i2]
                fExt = id + "Forward.fastq"
                rExt = id + "Reverse.fastq"
                writeOutFq(read, readR, fExt, rExt)
                continue
            else:
                belowThres += 1
                writeOutFq(read, readR)
                continue

        else:
            unmatchedInd += 0
            print(i2 + " Doesnts equal " + i6)
            writeOutFq(read, readR)
            continue

