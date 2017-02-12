#!/usr/bin/env python
'''
V1.0
Updated on 2013-11-26
'''

from File_Class import *
from General_Library import *
from Constant_Library import *

def identify_polyU(SEQ):
    sequence_length=len(SEQ)
    polyU_len=0
    for i in range(sequence_length-1,-1,-1):
        nt=SEQ[i]
        if nt=="T":
            polyU_len+=1
        else:
            break
    if polyU_len==0:
        trimmed_SEQ=SEQ
    else:
        right_border=polyU_len*(-1)
        trimmed_SEQ=SEQ[:right_border]
    return trimmed_SEQ,polyU_len

def identify_polyU_allowMismatch(SEQ,MAX_MISMATCH_COUNT=1,MAX_MISMATCH_RATIO=0.8):
    sequence_length=len(SEQ)
    polyU_len=0
    mismatch_count=0
    for i in range(sequence_length-1,-1,-1):
        nt=SEQ[i]
        if nt=="T":
            polyU_len+=1
        else:
            mismatch_count+=1
            mismatch_ratio=float(mismatch_count)/(mismatch_count+polyU_len)
            if mismatch_count>MAX_MISMATCH_COUNT:
                break

    return polyU_len

def identify_polyU_allowMismatch_seq(SEQ,MAX_MISMATCH_COUNT=1,MAX_MISMATCH_RATIO=0.8):
    ## this function still does not handle with mismatch section very well
    ## therefore, still needs to be improved for mismatch cases
    sequence_length=len(SEQ)
    polyU_len=0
    mismatch_count=0
    for i in range(sequence_length-1,-1,-1):
        nt=SEQ[i]
        if nt=="T":
            polyU_len+=1
        else:
            mismatch_count+=1
            mismatch_ratio=float(mismatch_count)/(mismatch_count+polyU_len)
            if mismatch_count>MAX_MISMATCH_COUNT:
                break


    if polyU_len>0:
        trimmed_SEQ=SEQ[:(-1)*(polyU_len+MAX_MISMATCH_COUNT)]
    else:
        trimmed_SEQ=SEQ
    return trimmed_SEQ,polyU_len

def find_ID_from_filename(filename):
    infile_name_list=filename.split("_")
    infile_name_length=len(infile_name_list)
    for index in xrange(infile_name_length):
        name_info=infile_name_list[index]
        if name_info[0]=="a" and (name_info[1] in digit_list):
            final_ref_ID=name_info+";"+infile_name_list[index+1]
    return final_ref_ID

def search_longest_T(seq,polyT_length_min=4):
    total_length=len(seq)
    longest_polyT_seq="T"
    current_seq=""
    for index in range(total_length):
        nt=seq[index]
        if nt=="T":
            current_seq=current_seq+"T"
        else:
            if current_seq!="":
                if len(current_seq)>len(longest_polyT_seq):
                    longest_polyT_seq=current_seq
                current_seq=""
    if len(current_seq)>len(longest_polyT_seq):
        longest_polyT_seq=current_seq
    max_length=len(longest_polyT_seq)
    max_coor=seq.find(longest_polyT_seq)
    max_end_coor=max_coor+max_length
    unpolyU_end_coor=max_coor-1
    return max_length,max_end_coor,unpolyU_end_coor

def search_longest_A(seq,polyA_length_min=4):
    total_length=len(seq)
    longest_polyA_seq="A"
    current_seq=""
    for index in range(total_length):
        nt=seq[index]
        if nt=="A":
            current_seq=current_seq+"A"
        else:
            if current_seq!="":
                if len(current_seq)>len(longest_polyA_seq):
                    longest_polyA_seq=current_seq
                current_seq=""
    if len(current_seq)>len(longest_polyA_seq):
        longest_polyA_seq=current_seq
    max_length=len(longest_polyA_seq)
    max_coor=seq.find(longest_polyA_seq)
    max_end_coor=max_coor+max_length
    unpolyU_end_coor=max_coor-1
    return max_length,max_end_coor,unpolyA_end_coor

def read_cluster_output(cls_output):
    cls_obj=GeneralFile_class(cls_output)
    cls_reader=cls_obj.reader_gen()
    row_count=0
    cluster_dict=dict()
    for row in cls_reader:
        row_count+=1
        if row[0][0]==">":
            if row_count==1:
                cluster_list=[]
            else:
                cluster_dict[cluster_id]=cluster_list
                cluster_list=[]
        else:
            #debug
            seq_info=row[1]
            cluster_id_list=seq_info.split(">")
            cluster_id_info=cluster_id_list[1]
            seq_id=""
            if cluster_id_info.count(".")==4:   
                for item in cluster_id_info:
                    if item!=".":
                        seq_id=seq_id+item
                    else:
                        break
            elif cluster_id_info.count(".")==4 and cluster_id_info.count("*")==1:
                sep_list = cluster_id_info.split(".")
                seq_id = sep_list[0]+"."+sep_list[1]

            elif cluster_id_info.count(".")==5 and cluster_id_info.count("*")==0:
                sep_list = cluster_id_info.split(".")
                seq_id = sep_list[0]+"."+sep_list[1]

            else:
                print cluster_id_info
                print "format unsupported"
                sys.exit(0)

            

            cluster_list.append(seq_id)

            if seq_info.count("*")==1:    
                cluster_id=seq_id
    ##out the last one
    cluster_dict[cluster_id]=cluster_list
    return cluster_dict