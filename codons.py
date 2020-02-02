#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:04:15 2020

@author: alexander g. lucaci


The idea for this is that:
    
    Given a protein sequence and a transcript sequence
    I find the codons by stepping over the transcript sequence until the translated sequence matches the protein sequence
    that way, I have only the codons and not the additional sequences from the transcript
    (Which may be useful later)
    I will also create two output files
        One with the STOP codon stripped (this makes it hyphy compatible.)
        One with the STOP codons (may be useful later, codon bias?)
"""

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
from Bio.Alphabet import generic_rna

# =============================================================================
# Declares
# =============================================================================
PROTEIN = "TP53_refseq_protein.fasta"
TRANSCRIPTS = "TP53_refseq_transcript.fasta"
OUTPUT = "TP53_refseq_CODONS.fasta"


# =============================================================================
# Helper functions
# =============================================================================
#turn into class

def Process(PROTEIN_DESC, PROTEIN_SEQ, TRANSCRIPT_DESC, TRANSCRIPT_SEQ):
    #print([PROTEIN_DESC], [TRANSCRIPT_DESC])
    print("Protein seq length (AA):", len(PROTEIN_SEQ))
    print(PROTEIN_DESC, PROTEIN_SEQ)
    
    
    TRIMMED_TRANSCRIPT_SEQ = ""
    print()
    
    #Loop over TRANSCRIPT_SEQ
    start = 0
    NT_SEQ_LENGTH = len(PROTEIN_SEQ) * 3
    while start < len(str(TRANSCRIPT_SEQ)):
        
        coding_dna = TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH].translate()
        
        
        #if start == 202:
        #    print(start, coding_dna, "\n", len(coding_dna))
        
        if coding_dna == str(PROTEIN_SEQ):
            print("\n#### FOUND", coding_dna)
            #print("#### CODONS", TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH + 3]) # has stop codon
            print("\n#### CODONS", TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH]) # NO stop codon
            break
        
        start += 1
        if start == 301: break
        #end if
    #end while
    
    return TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH]

# =============================================================================
# Main subroutine.
# =============================================================================



def main(): # Really to verify things.
    global PROTEIN, TRANSCRIPTS
    print("\tTRANSCRIPT INPUT FILE:", TRANSCRIPTS)
    print("\tPROTEIN INPUT FILE:", PROTEIN)
    print()
    
    protein_list = []
    transcript_list = []
    
    with open(TRANSCRIPTS, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            count +=1
            transcript_list.append(record.description)
        print("\tTranscripts:", count)    
    handle.close()
    
    with open(PROTEIN, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            count +=1
            protein_list.append(record.description)
        print("\tProteins:", count)
    handle.close()
    
    
    #for n, item in enumerate(transcript_list):
    #    print(item, [protein_list[n]])
    
# =============================================================================
# Main
# =============================================================================
#Verify files exist
print("# =============================================================================")
print("# Processing... ")
main()
print("# =============================================================================")
#Looks like species all match up in transcript and protein fasta.
#This is exceptional, will need to look for species name (from protein desc.) in transcript desc.

with open(OUTPUT, "w") as fh:
    fh.write("")
fh.close()

with open(PROTEIN, "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, "fasta")):
        #if n == 1: break
        print("\n" + str(n))
        protein_id = record.id
        protein_desc = record.description
        protein_seq = record.seq
        
        
        with open(TRANSCRIPTS, "r") as transcript_handle:
            for m, transcript_record in enumerate(SeqIO.parse(transcript_handle, "fasta")):
                if m == n:
                    transcript_id = transcript_record.id
                    transcript_desc = transcript_record.description
                    transcript_seq = transcript_record.seq
                #end if
            #end inner for
        transcript_handle.close()
        #end inner with
        
        
        #Process
        codons = Process(protein_desc, protein_seq, transcript_desc, transcript_seq)
        
        #Print out transcript desc and TRIMMED codons transcript.
        with open(OUTPUT, "a") as fh:
            fh.write(">" + transcript_desc + "\n" + str(codons) + "\n")
        fh.close()
        
    #end outer for
prot_handle.close()


#end outer with



                    
                
                
                
                
                
                
        


# =============================================================================
# End of file    
# =============================================================================