## Demultiplexing Assignment
Code developed for paired end reads in which reads and indexes are in seperate files. The following scripts
will output write statistics and a graph of relative sample abundances for the demultiplexing. 

__Usage__ 
DeMultiAlgoArg.py [options]

__Input flags__

_Required_

-r1 [readfile1] 

-i1 [indexfile1] 

-r2 [readfile2] 

-i2 [indexFile2] 

_Optional_

-t [mean quality score threshold] _Default:_ 25

-b [Barcodes mapping file] _Default:_ indexes.txt


### Statistics defined

__Records Read:__ Total number of reads processed in each file. All four files should have the same amount of lines, if not, the script will only read until your shortest file.

__bad Reads:__ N found in barcode of either of the 2 index reads

__Below Threshold:__ Average quality scores of each index read does not surpass threshold. No N's in barcode though

__unmatched indices:__ Indexes don't match. Passes threshold and N tests stated above.

__Not found in indexes file:__ Index of forward read not found in indexes.txt input file. Passes all tests above. 

__Matched Reads:__ Passes all above test. Read written to respective sample with forward and reverse written seperately. Forward index place inline with read headers for forward and reverse reads.