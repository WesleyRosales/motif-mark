#!/usr/bin/env python

import cairo
import math
import argparse

def get_args():
    my_parser = argparse.ArgumentParser(prog="Motif Mark",
                description='Outputs visual representation of input DNA and input motifs found within DNA seqeunce.')
    
    my_parser.add_argument('-f',
                        '--fasta',
                        action='store',
                        help='Enter DNA fasta file',
                        required=True)

    my_parser.add_argument('-m',
                        '--motif',
                        action='store',
                        help='Enter IUPAC notation motif file',
                        required=True)

    return my_parser.parse_args()

args = get_args()


# function to translate IUPAC nucleotide notation to DNA/RNA sequences
# input: IUPAC notation sequence
# output: list of all DNA/RNA sequences corresponding to that notation
def iupac_convert(seq: str) -> list:

    seq = seq.upper()

    conversion_dict = {
        'A':['A'],
        'C':['C'],
        'G':['G'],
        'T':['T'],
        'U':['U'],
        'W':['A','T'],
        'S':['C','G'],
        'M':['A','C'],
        'K':['G','T'],
        'R':['A','G'],
        'Y':['C','T'],
        'B':['C','G','T'],
        'D':['A','G','T'],
        'H':['A','C','T'],
        'V':['A','C','T'],
        'N':['A','C','T'],
        '-':[]
    }

    translated_list = []
    for letter in seq:
        translated_list.append(conversion_dict[letter])

    output_list = []

    for lst in translated_list:
        new_list = []
        if len(output_list) == 0: #creates first letter of motifs in list
            for k in range(len(lst)):
                output_list.append(lst[k])
            continue
        for j in range(len(lst)):
            working_list = output_list
            for i in range(len(working_list)):
                new_string = working_list[i] + lst[j] #adds the element to each string in the output list
                new_list.append(new_string)

        output_list = new_list[:]

    return output_list 

# finds the start positions of introns of an input sequence
# input must be formatted so that the introns are lowercase and exons are uppercase
def gene_architecture_positions(dna: str) -> list:
     
    count = -1 #create count that matches index-0 for future use in defining position of intron/exon in sequence
    positions = [[],[]] # nested lists where first list = intron start positions and second list = exon start positions
    for ch in dna:
        count += 1
        
        # see if beginning of sequence is intron or exon
        if count == 0:
            if ch.islower() == True:
                currently_lower = True
                positions[0].append(count)
            else:
                currently_lower = False
                position[1].append(count)

        # use state of currently_lower to find where the state of the cases changes and add position accordingly
        if ch.islower() == True and currently_lower == False:
            currently_lower = True
            positions[0].append(count)
        
        if ch.islower() == False and currently_lower == True:
            currently_lower = False
            positions[1].append(count)

    return positions

# save all motifs within a list


# motif object to store motif as input and the specific sequences corresponding to that motif
class mtf:
    def __init__(self, iupac):
        
        ## Data ##
        self.iupac = iupac
        self.motifs = iupac_convert(iupac)



# gene object that stores:
    # name of gene
    # sequence length (base-0)
    # [[intron start positions],[exon start positions]]
    # motif start positions including which motif
class gene:
    def __init__(self, gene_name, sequence):
        
        ## Data ##
        self.gene = gene_name
        self.sequence = len(sequence) - 1
        self.intron_exon_start = gene_architecture_positions(sequence)
        

        ## Methods ##



#figure object that holds information relevent for cairo
    #gene name same as that of the gene class
    #intron and exon start positions to make the base of the figure
    #number of unique motifs found in gene class to determine number of tracks to put in figure
class figure:
    def __init__(self, gene_name):

        ## Data ##
        self.name = gene_name
        self.intron_exon_positions = gene_name.intron_exon_start

        ## Methods ##
        






print(iupac_convert("y"))

print(gene_architecture_positions("ggggGGGGgggg"))