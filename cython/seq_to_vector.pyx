#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:27:37 2017

@author: hilmar
"""

import cython
import numpy as np
cimport numpy as np

cdef unsigned int[256] nuc_char_to_code;
nuc_chars = "ACGTNUWSMKRYBDHV"
nuc_codes = "1234545555555555"

def set_translation():
    global nuc_char_to_code
    global nuc_codes
    global nuc_chars
    cdef int code, r, c
    
    for i in range(0,len(nuc_char_to_code)):
        nuc_char_to_code[i] = 0
    
    for r in range(0,len(nuc_chars)):
        c = ord(nuc_chars[r])
        code = int(nuc_codes[r])
        if c < 0 or c > 255:
            pass
        nuc_char_to_code[c] = int(code)
        

# FIXME: seq_startspos is not required here - remove on refactoring
def seq2num_vec(char* sequence, qual, long int seq_startpos, cigar_tuples, int alignment_len_default):
    """
    Convert read alignment to numeric vector with length equal to alignment length (plus soft clipped bases)
    """
    
    cdef unsigned int cur_alignment_vec_len = alignment_len_default
    cdef unsigned int cig_pos = 0 # this is the position after the last inserted cigar operation relativ to the leftmost alignment start on the ref seq
    cdef unsigned int q_pos = 0 # this is the corresponding position in the query sequence (i.e. the read sequence)
    cdef bint leftmost_match_seen = False
    cdef int start_position_offset = 0 # if the alignment starts with a soft clip then the start position of the read refers to the first really aligned base
    cdef unsigned int ll, i, ci, qi
    cdef int op_code
        
    #print(sequence)
    
    seq_vec = np.zeros(cur_alignment_vec_len, dtype=int)
    qual_vec = np.zeros(cur_alignment_vec_len, dtype=int)
    ins_list = {"I": {},"D":{}} # list of indel operations in the alignment
    #print cigar_tuples
    for cig_op in cigar_tuples:
        ll = cig_op[1]
        #print cig_op
        # check if the alignment vector is still long enough to store the remaining cigar op. If not, resize it
        if ll+cig_pos > cur_alignment_vec_len:
            seq_vec = np.append(seq_vec, np.zeros(cur_alignment_vec_len+ll, dtype=int))
            qual_vec = np.append(qual_vec, np.zeros(cur_alignment_vec_len+ll, dtype=int))                
            cur_alignment_vec_len += (cur_alignment_vec_len+ll)
        
        op_code = cig_op[0]
        #print op_code
        #print "ll: %d\tqpos: %d\tcig_pos: %d" % (ll, q_pos, cig_pos)
        if op_code in (0,7,8): # match
            #seq_vec[cig_pos:(cig_pos+ll)] = [int(e) for e in string.translate(sequence[q_pos:(q_pos+ll)],self.trans_tab)]
            #seq_vec[cig_pos:(cig_pos+ll)] = [nuc_char_to_code[ord(e)] for e in sequence[q_pos:(q_pos+ll)]]
            #qual_vec[cig_pos:(cig_pos+ll)] = qual[q_pos:(q_pos+ll)]
            for i in range(ll):
                qi = <unsigned int>(i+q_pos)
                ci = <unsigned int>(i+cig_pos)
                seq_vec[ci] = nuc_char_to_code[<unsigned int> sequence[qi]]
                qual_vec[ci] = qual[qi]

            cig_pos += ll # move the pointer in the result vector by cigar length
            q_pos += ll # move the pointer in the query (read) sequence by cigar length
            leftmost_match_seen = True
        elif op_code == 4: # soft clipping
            for i in range(ll):
                seq_vec[i+cig_pos] = -1
                qual_vec[i+cig_pos] = qual[i+q_pos]
            
            #seq_vec[cig_pos:(cig_pos+ll)] = -1
            #qual_vec[cig_pos:(cig_pos+ll)] = qual[q_pos:(q_pos+ll)]
            cig_pos += ll
            q_pos += ll
            if not leftmost_match_seen:
                start_position_offset += ll
        elif op_code == 2: # deletion
            for i in range(ll):
                seq_vec[i+cig_pos] = 6
                #qual_vec[i+cig_pos] = qual[i+q_pos]
                
            #seq_vec[cig_pos:(cig_pos+ll)] = 6
            ins_list["D"][q_pos] = ( ll, "" ) 
            cig_pos += ll
            leftmost_match_seen = True # It is not clear from the SAM spec if we could have xSyD and how that would affect start position. Lets assume that Del counts as alignment
        elif op_code == 1: # insertion
            #print ">%s<" % sequence[q_pos:(q_pos+ll)]
            ins_list["I"][q_pos] = ( ll, sequence[q_pos:(q_pos+ll)] ) 
            q_pos += ll
            leftmost_match_seen = True # It is not clear from the SAM spec if we could have xSyI and how that would affect start position. Lets assume that Ins counts as alignment
        elif op_code in (5,6): # hard clipping, padding
            continue
    return (seq_vec, qual_vec, ins_list, cig_pos, start_position_offset)
