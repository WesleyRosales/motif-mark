#!/usr/bin/env python

import cairo
import math
import argparse
import re

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
        'T':['T','U'],
        'U':['U','T'],
        'W':['A','T','U'],
        'S':['C','G'],
        'M':['A','C'],
        'K':['G','T','U'],
        'R':['A','G'],
        'Y':['C','T','U'],
        'B':['C','G','T','U'],
        'D':['A','G','T','U'],
        'H':['A','C','T','U'],
        'V':['A','C','T','U'],
        'N':['A','C','T','U'],
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

# finds the start and stop positions of exons
# input must be formatted so that the introns are lowercase and exons are uppercase
def gene_architecture_positions(dna: str) -> list:
     
    count = -1 #create count that matches index-0 for future use in defining position of intron/exon in sequence
    positions = [[],[]] # nested lists where first list = intron start positions and second list = exon start positions
    currently_lower = True
    for ch in dna:
        count += 1
        
       # use state of currently_lower to find where the state of the cases changes and add position accordingly
        if ch.islower() == False and currently_lower == True:
            currently_lower = False
            positions[0].append(count)        
        
        if ch.islower() == True and currently_lower == False:
            currently_lower = True
            positions[1].append(count-1) # count - 1 is true end of exon
        
        if count == len(dna) - 1 and ch.islower == False:
            positions[1].append(count)
    return positions


# motif object to store motif as input and the specific sequences corresponding to that motif
class Mtf:
    def __init__(self, iupac):
        ## Data ##
        self.iupac = iupac
        self.motifs = iupac_convert(iupac)



# gene object that stores:
    # name of gene
    # sequence length
    # [[intron start positions],[exon start positions]]
    # motif start positions including which motif
class Gene:
    def __init__(self, gene_name, sequence):
        ## Data ##
        self.gene = gene_name
        self.sequence = sequence
        self.exon_start_stop = gene_architecture_positions(sequence)
        self.sequence_length = len(sequence)
        
    ## Methods ##
    # uses regex to find all unique start locations of an input motif
    def motif_locations(self, mtf_object):
        start_stop = [[],[]]
        for motif in mtf_object.motifs:
            x = re.finditer(motif, self.sequence, re.IGNORECASE)
            if x != None:
                for match in x:
                    start_stop[0].append(match.start())
                    start_stop[1].append(match.start()+len(motif)-1)
        
        return start_stop

            


#figure object that holds information relevent for cairo
    #gene name same as that of the gene class
    #intron and exon start positions to make the base of the figure
    #number of unique motifs found in gene class to determine number of tracks to put in figure
class Figure:
    def __init__(self, gene_name, width, height):

        ## Data ##
        self.name = gene_name
        self.exon_start_stop = gene_holder[gene_name].exon_start_stop
        self.height = height # indicates where figure starts on overall surface (top-left corner)
        self.width = width # consistent for all figures

    ## Methods ##
    
    #create the base for each figure including: intron, exon, figure separation line, header
    def base(self):
        ctx.set_source_rgb(0, 0, 0)
        # header in top-left corner of figures
        ctx.select_font_face("Arial",
                    cairo.FONT_SLANT_NORMAL,
                    cairo.FONT_WEIGHT_BOLD)
        ctx.set_font_size(15)
        ctx.move_to(0, self.height)
        ctx.show_text(self.name)
        # intron line
        ctx.set_line_width(5)
        ctx.move_to(0, self.height + 250)
        ctx.line_to(self.width, self.height + 250)
        ctx.stroke()
        #bottom border of figure
        ctx.set_line_width(2)
        ctx.move_to(0, self.height + 280)
        ctx.line_to(self.width, self.height + 280)
        ctx.stroke()
        # exons
        ctx.set_line_width(20)
        for i in range(len(self.exon_start_stop[0])):
            ctx.move_to(self.exon_start_stop[0][i], self.height + 250)
            ctx.line_to(self.exon_start_stop[1][i], self.height + 250)
            ctx.stroke()
    
    def find_motifs(self):
        color_count = 0
        figure_scale_factor = self.width/gene_holder[self.name].sequence_length
        ctx.set_line_width(15)
        ctx.set_source_rgb(color_palette[3][0], color_palette[3][1], color_palette[3][2])

        for motif in motif_holder.keys():
            start_stop = gene_holder[self.name].motif_locations(motif_holder[motif])
            print(start_stop)
            ctx.set_source_rgba(color_palette[color_count][0], color_palette[color_count][1], color_palette[color_count][2])
            color_count += 1
            for i in range(len(start_stop[0])):
                ctx.move_to(start_stop[0][i]*figure_scale_factor, self.height + 250 - (20*(color_count)))
                ctx.line_to(start_stop[1][i]*figure_scale_factor, self.height + 250 - (20*(color_count)))
                ctx.stroke()
    





# save all motifs within a list
motif_file = open(args.motif, "r")

motif_list = motif_file.read().splitlines()

motif_file.close()

#dictionary containing motif objects
motif_holder = {name: Mtf(name) for name in motif_list}


#gene_holder will be a dictionary of with the name of the sequence as the key and the gene object as the value
gene_holder = {}
new_string = ""

seq_file = open(args.fasta, "r")

# line-by-line, look for header and sequence. sequence lines will be stripped on \n and concatenated to a new string.

for line in seq_file:
    line = line.strip()
    if line not in gene_holder.keys() and line.startswith(">"):
        if len(new_string) != 0:
            gene_holder[header] = Gene(header, new_string)
            new_string = ""
        header = line
        continue
    new_string += line

gene_holder[header] = Gene(header, new_string) # add the last entry after the last line is processed

seq_file.close()


# width of figure stays constant; height = (# of motifs * 20) + (number of sequences * 300)
heading = (len(motif_list) + 2)*20
body = len(gene_holder.keys())*300
WIDTH, HEIGHT = 1000, heading + body

# dictionary to hold figure objects
figure_holder = {}
fig_count = 0
for key in gene_holder.keys():
    figure_holder[key] = Figure(key, WIDTH, heading + fig_count*300)
    fig_count += 1


#cairo.ImageSurface is used for pngs
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
ctx = cairo.Context(surface)

# BACKGROUND
ctx.set_source_rgb(0.96, 0.96, 0.76) #beige
ctx.rectangle(0, 0, WIDTH, HEIGHT)  # Rectangle(x0, y0, x1, y1) of entire width and height
ctx.fill() # fill the rectangle

# COLOR PALETTE
#length 10 list of colors to use for motifs
# NOTE: custom palette, may not be color-blind friendly. should review before submission
color_palette = [[221/255, 106/255, 173/255],
                [0, 105/255, 137,255],
                [124/255, 11/255, 43/255],
                [127/255, 41/255, 130/255],
                [117/255, 142/255, 79/255],
                [234/255, 191/255, 203/255],
                [226/255, 109/255, 90/255],
                [17/255, 157/255, 164/255],
                [204/255, 0/255, 122/255],
                [160/255, 193/255, 209/255]]

# MOTIF LEGEND
ctx.set_source_rgb(0, 0, 0) #set the color of title to black

ctx.set_font_size(15)
#may want to change font before submitting
ctx.select_font_face("Arial",
                    cairo.FONT_SLANT_NORMAL,
                    cairo.FONT_WEIGHT_BOLD)

#move the text to the top-right corner of the surface
ctx.move_to(WIDTH*(8/10), 15)
ctx.show_text("MOTIF") 

#cycle through motif_list and print each one below the previous
for j in range(len(motif_list)):
    ctx.set_source_rgba(color_palette[j][0], color_palette[j][1], color_palette[j][2])
    track = (j+2)*20
    ctx.move_to(WIDTH*(8/10), track)
    ctx.show_text(motif_list[j])

#calling the methods in the figure objects to create the figure
for key in figure_holder.keys():
    print(key)
    figure_holder[key].base()
    figure_holder[key].find_motifs()

surface.write_to_png("TEST_DONOTKEEP.png")  # Output to PNG


#NOTE: plan going forward is to find a way to create all the objects and use them to create
# the figure.
# 4. alter color scheme to make each motif more identifiable.



