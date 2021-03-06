# Motif Mark

Proteins may interact with DNA by recognizing stretches of DNA called "motifs". These motifs can be highly specific (such as ATGCCA) or more loose (such as YGRAN). The human eye may be able to recognize these motifs in a short DNA sequence however, it can become overwhelming when there are tens of motifs to search for in a sequence measured in kilobases. 
<br>
<br>
The purpose of Motif Mark is to create a visual representation of DNA sequences with correct anatomy (differentiate between exons and introns). Along that representation, motifs will be identified and labelled. Using this representation instead of raw DNA sequences will allow the user to visually inspect their DNA input for their motifs and where they fall within the input DNA.
<br>
<br>
Motif Mark takes in two inputs: a fasta file (designated as -f), and a text file containing a single IUPAC-notation motif per line. The data folder contains example data. It will then output a single png image containing motif figures for each sequence within the fasta file at the location motif-mark-oop.py is run.
<br>
<br>
Example command line for code (from withing motif-mark folder)
```
python motif-mark-oop.py -f <FASTA FILE> -m <MOTIF FILE>
```