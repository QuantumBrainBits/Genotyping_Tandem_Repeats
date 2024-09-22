# packages

import re
from tqdm import tqdm
from ast import literal_eval
import os
import gzip
from tqdm import tqdm
from collections import Counter
import os
import numpy as np
import itertools
import re
import numpy as np
import pysam
from pysam import VariantFile
import json
from statistics import mean
import statistics
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import time
import tracemalloc



class realign:

    # default running def and initilization of variables.
    def __init__(self, duplicated_read, read_reference_start, read_reference_end, query_seq, query_id, qualities_, bases_, cigar_len, cigar_type, end_pos_, pos_):
        
        self.duplicated_read       = duplicated_read
        self.read_reference_start  = read_reference_start
        self.read_reference_end    = read_reference_end
        self.query_seq             = query_seq
        self.query_id              = query_id
        self.qualities_            = qualities_
        self.bases_                = bases_
        self.cigar_len             = cigar_len
        self.cigar_type            = cigar_type
        self.end_pos_              = end_pos_
        self.pos_                  = pos_
        # self.recursion             = recursion
        # self.Trimmed_low_quals_list= Trimmed_low_quals_list
    
    # 
    def TrimAlignment(self, min_read_start, max_read_stop, min_base_qual, Trim_qual_bases):

            
        ltrim     = 0
        start_pos = self.pos_

        # As we use the recursion for both "Ends trimming" & "Far away STRs trimming".
        if Trim_qual_bases == True: # when trimming "low qual bases"
            i = 1
        else: i = 0 # this is when "Trim bases which are far away"
        
        # left to right trimming.        
        while((start_pos < min_read_start) and (len(self.cigar_type) > i)): #and len(cigar_ops_) > 0) ---note: cigar len is being checked & i assume all cigar lens will be greater than 0, until it is unmapped.
            
            #check if we should stop trimming b/c the quality score is above the threshold.
            qual_above_thresh = False

            # Trim the 40bp either sides of the STR.
            if Trim_qual_bases == True:    
                if self.cigar_type[0] in ['M', '=', 'X', 'I', 'S']: # do notin if 'D', 'H', cigar_ops[ltrim] was used 
                         qual_above_thresh = (self.qualities_[ltrim] > min_base_qual)
                         
            # here we check if qual is above threshold
            if qual_above_thresh == True:
                break
             
            # here we increament the start and ltrim positions. # do notin if 'H'
            if self.cigar_type[0] in ['M', '=', 'X']:
                ltrim += 1
                start_pos += 1
            elif self.cigar_type[0] in ['D']: start_pos += 1
            elif self.cigar_type[0] in ['I', 'S']: ltrim += 1

            # removing the bases using CIGAR string.
            if self.cigar_len[0] == 1:
                self.cigar_len.remove(cigar_len[0])
                self.cigar_type.remove(cigar_type[0])
            else: self.cigar_len[0] -= 1

            # removing the base from the read sequence, ( only I, S, M can be removed, D is not present in the read.)
            

        # Right to left trimming.
        rtrim = 0
        qual_string_len = int(len(self.qualities_)) - 1
        rtrim_index     = -rtrim-1
        end_pos         = self.end_pos_

        while((end_pos > max_read_stop) and (len(self.cigar_type) > i)):


            self.qual_above_thresh = False

            # Trim the 40bp either sides of the STR.
            if Trim_qual_bases == True:    
                #
                if  self.cigar_type[0] in ['M', '=', 'X', 'I', 'S']:
                    self.qual_above_thresh = (self.qualities_[rtrim] > min_base_qual)

            
            # here we check if qual is above threshold
            if self.qual_above_thresh == True: break

            # here we increament the start and rtrim positions.
            if self.cigar_type[0] in ['M', '=', 'X']:
                rtrim -= 1
                end_pos -= 1
            elif self.cigar_type[0] in ['D']: end_pos -= 1
            elif self.cigar_type[0] in ['I', 'S']: rtrim -= 1

            # removing the bases using CIGAR string.
            if self.cigar_len[-1] == 1:
                self.cigar_len.remove(self.cigar_len[-1])
            else: self.cigar_len[-1] -= 1
                
                        

        # 
        assert(ltrim + rtrim <= len(self.bases_)), "Checking if all bases are trimmed"
        self.bases_     = self.bases_[ltrim : len(self.bases_)+rtrim]
        self.qualities_ = self.qualities_
        length_         = ltrim - rtrim
        self.pos_       = start_pos
        self.end_pos_   = end_pos


        return self.bases_, start_pos, end_pos, self.cigar_len, self.cigar_type

    
    # Function pases the min base qual and returns list with deleted low base qualities.
    def TrimLowQualityEnds(self, min_base_qual):

        return realign.TrimAlignment(self, self.end_pos_+1, self.pos_-1, min_base_qual, True)
