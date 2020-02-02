#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 12:40:07 2020

@author: alexander g. lucaci


The idea for this is to filter out "PREDICTED" or "LOW QUALITY Sequences"
"""
# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO

# =============================================================================
# Declares
# =============================================================================
PROTEIN = "TP53_refseq_protein.fasta"
TRANSCRIPTS = "TP53_refseq_transcript.fasta"



# =============================================================================
# Main
# =============================================================================



with open(TRANSCRIPTS, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if "PREDICTED" not in record.description:
            #print(record, "\n")
            #print(record.description)
            print(">" + str(record.description), record.seq)
        
# =============================================================================
# End of file    
# =============================================================================
