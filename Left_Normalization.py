# realign the reads with left normalization concept. (AlignmentOps)
'''
Realign read to reference region using left alignment variant.
Store the new alignment info using the provided alignment reference.
Converts the bases to their upper case variants.
'''
# fetch part of the sequence "sequence 1" from position 0 to 3
from pysam import FastaFile

# Do you know what is this cell doing ? : 1.

def Generate_realign_seq(chrom, ref_start, ref_end):

    fasta = "../../../../malini/repeats/reference/Homo_sapiens_assembly38.fasta"
    sequences_object   = FastaFile(fasta)
    # chrom_size         = len(sequences_object.fetch(chrom)
    reference_sequence = sequences_object.fetch(chrom, ref_start, ref_end)

    # return
    return reference_sequence


# function which creates the candidate seqs 
# def Gen_candidate_seqs():

    
    