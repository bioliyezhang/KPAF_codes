#!/usr/bin/env python
'''
V2015-05-02
'''
import os
from File_Class import *
from General_Library import *
from Constant_Library import *
import numpy

'''
def updates_plans(documents):
    #
    #
    return Null

'''

'''
def updates_notes(documents):
    #updates:2015-05-02
    #1. add new function to extract count info
    return 0

'''

'''
def fast_only_retrieve_bases(read,cigar,read_start_coor,search_coor):
        
        
    ## this will provide the two border of filtered range
    #filtered_out=False
    #SNV_base=""
    read_length=len(read)
    #print "READ Length", read_length
    #print "read length",read_length
    
    #final_data_list=[]
    temp_base=""
    mapping_pattern_list=[]
    border_list=[]
    read_left_set=False
    read_right_set=False
    search_coor_found=False
    read_right_coor=0
    #mapping_count_list=[]
    #M_count=0
    for character in cigar:
        if character in digit_list:
            temp_base=temp_base+character
        else:
            #fragment_length=int(temp_base)
            
            fragment_data=temp_base+character
            temp_base=""
            #final_data_list+=fragment_length*[character]
            mapping_pattern_list.append(fragment_data)
            #mapping_count_list.append(fragment_length)
                
    
    ##Process the data
    current_read_index=-1
    current_coor_index=read_start_coor-1
    previous_pattern=""
    
    for data in mapping_pattern_list:
        data_length=int(data[:-1])
        data_type=data[-1]
        if data_type=="M":
            current_coor_index+=data_length
            current_read_index+=data_length
        elif data_type=="N":
            border_list.append(current_coor_index)
            current_coor_index+=data_length
            border_list.append(current_coor_index)
        elif data_type=="I":
            current_read_index+=data_length
        elif data_type=="D":
            current_coor_index+=data_length
        else:
            print "incorrect data format,check input"
            sys.exit(0)
        
            
        if current_read_index>=5 and read_left_set==False:
            read_left_coor=current_coor_index-(current_read_index-5)
            read_left_set=True
        
        if current_read_index>=(read_length-6) and read_right_set==False:
            read_right_coor=current_coor_index-(current_read_index-(read_length-6))
            read_right_set=True
        
        if current_coor_index>=search_coor and search_coor_found==False:
            #search_coor_mark=data_type
            if data_type=="N":
                SNV_base="N"
            elif data_type=="N":
                SNV_base="D"
            else:
                current_read_index=current_read_index-(current_coor_index-search_coor)
                SNV_base=read[current_read_index]
            search_coor_found=True
            break
    
    near_splice_junction=False  
    for border_coor in border_list:
        border_distance=abs(border_coor-search_coor)
        if border_coor<=10:
            near_splice_junction=True
            break
    
    #print SNV_base
    return SNV_base 
'''

def read_cigar(cigar_data):
    #print cigar_data
    result = []
    match_count = 0
    cigar_length=len(cigar_data)
    tmp=""
    for index in range(cigar_length):
        item = cigar_data[index]
        if item in CIGAR_LIST:
            output = tmp+item
            result.append(output)
            if item=="M":
                match_count+=int(tmp)
            tmp=""
        elif item in DIGIT_LIST:
            tmp=tmp+item
        else:
            print "encounter non-canonical CIGAR string"
            print "the non-canonical CIGAR charater is", item
            print "the non-canonical CIGAR string is", cigar_data
            sys.exit(0)
    return result,match_count

def read_cigar_withoutID(cigar_data):
    #print cigar_data
    result = []
    match_count = 0
    cigar_length=len(cigar_data)
    tmp=""
    for index in range(cigar_length):
        item = cigar_data[index]
        if item in CIGAR_LIST:
            output = tmp+item
            if item=="I" or item=="D":
                pass
            else:
                result.append(output)
            if item=="M":
                match_count+=int(tmp)
            tmp=""
        elif item in DIGIT_LIST:
            tmp=tmp+item
        else:
            print "encounter non-canonical CIGAR string"
            print "the non-canonical CIGAR charater is", item
            print "the non-canonical CIGAR string is", cigar_data
            sys.exit(0)
    return result,match_count

def read_cigar_with_length(cigar_data):
    
    result = []
    match_count = 0
    total_length = 0
    cigar_length=len(cigar_data)
    tmp=""
    for index in range(cigar_length):
        item = cigar_data[index]
        if item in CIGAR_LIST:
            output = tmp+item
            result.append(output)
            if item=="M":
                match_count+=int(tmp)
            if item!="D":
                total_length+=int(tmp)
            tmp=""
        elif item in DIGIT_LIST:
            tmp=tmp+item
        else:
            print "encounter non-canonical CIGAR string"
            print "the non-canonical CIGAR charater is", item
            print "the non-canonical CIGAR string is", cigar_data
            sys.exit(0)

    ## result is a list of cigar string with length as a item in the list
    ## match_count is the total matched length 
    ## total length is the total length of the fragment 
    return result,match_count,total_length

def read_cigar_coor_offset(cigar_data):
    #print cigar_data
    result = []
    match_count = 0
    total_length = 0
    cigar_length=len(cigar_data)
    tmp=""
    for index in range(cigar_length):
        item = cigar_data[index]
        if item in CIGAR_LIST:
            output = tmp+item
            result.append(output)
            if item=="M":
                match_count+=int(tmp)
            if item=="D":
                match_count+=int(tmp)
            if item!="I":
                total_length+=int(tmp)
            tmp=""
        elif item in DIGIT_LIST:
            tmp=tmp+item
        else:
            print "encounter non-canonical CIGAR string"
            print "the non-canonical CIGAR charater is", item
            print "the non-canonical CIGAR string is", cigar_data
            sys.exit(0)
    return total_length,match_count

def read_cigar_coor_offset_withlist(cigar_data):
    #print cigar_data
    result = []
    match_count = 0
    total_length = 0
    cigar_length=len(cigar_data)
    tmp=""
    for index in range(cigar_length):
        item = cigar_data[index]
        if item in CIGAR_LIST:
            output = tmp+item
            result.append(output)
            if item=="M":
                match_count+=int(tmp)
            if item=="D":
                match_count+=int(tmp)
            if item!="D":
                total_length+=int(tmp)
            tmp=""
        elif item in DIGIT_LIST:
            tmp=tmp+item
        else:
            print "encounter non-canonical CIGAR string"
            print "the non-canonical CIGAR charater is", item
            print "the non-canonical CIGAR string is", cigar_data
            sys.exit(0)
    return result,total_length,match_count


def reverse_cigar(cigar_data):
    '''
    version: 2015-03-30
    
    '''
    reversed_cigar =""
    cigar_length=len(cigar_data)
    tmp=""
    previous_cigar=""
    for index in range(cigar_length-1,-1,-1):
        item = cigar_data[index]
        if item in CIGAR_LIST:
            if previous_cigar!="":
                output = tmp+previous_cigar
                reversed_cigar=reversed_cigar+output
                tmp=""
            previous_cigar=item

        elif item in DIGIT_LIST:
            tmp=item+tmp
        else:
            print "encounter non-canonical CIGAR string"
            print "the non-canonical CIGAR charater is", item
            print "the non-canonical CIGAR string is", cigar_data
            sys.exit(0)
    
    output = tmp+previous_cigar
    reversed_cigar=reversed_cigar+output
    return reversed_cigar


def full_filter_retrieve_bases(read,cigar,read_start_coor,search_coor,read_quality,READEND_CUTOFF=6,SPLICEJUNCTION_CUTOFF=10):
    
    ##update on the deletion criteria: 2013.11.11
    
    ##[SectionI]
    ##Variable Initialization   
    read_length=len(read)
    filter_out=False    ##filter_out indicates whether this read confidently support the SNV calls
    temp_base=""        ##temp variable to store the cigar patterns 
    mapping_pattern_list=[]   ##This list stores the separated cigar patterns for READ alignments
    border_list=[]      ## This list stores the splice-junction coordinates 
    read_left_set=False ## this record whether left end read has been determined
    read_right_set=False    ## this record whether right end read has been determined
    search_coor_found=False ## this record whether SNV coordinate has been located
        
    ##[SectionII]
    ##Process the cigar string and separate it into Lists
    for character in cigar:
        if character in digit_list:   ## test whether the string is still a digit 0-9
            temp_base=temp_base+character     ##compile number 
        else:
            fragment_data=temp_base+character   ##other wise create cigar item
            temp_base=""
            mapping_pattern_list.append(fragment_data)
                
    ##[Section III]
    ##locate the coordinate bases
    current_read_index=-1           ## record the read location
    current_coor_index=read_start_coor-1    ## record the genome coordinate location
    SNV_base="N"
    
    for data in mapping_pattern_list:
        data_length=int(data[:-1])  ## get the cigar length 
        data_type=data[-1]      ## get the cigar type
        
        ##update the coordinate and read location
        if data_type=="M":      ## M means mapped
            current_coor_index+=data_length     ##update genome coor
            current_read_index+=data_length     ##update read coor
        elif data_type=="N":        ## N means splicing (jumpped length)
            border_list.append(current_coor_index)
            current_coor_index+=data_length
            border_list.append(current_coor_index+1)
        elif data_type=="I":        ## I means insertion
            current_read_index+=data_length
        elif data_type=="D":        ## D means deletion
            current_coor_index+=data_length
        elif data_type=="S":
            current_read_index+=data_length
        else:
            print "incorrect data format,check input"
            sys.exit(0)
        
        ## Setup left read border   
        if current_read_index>=(READEND_CUTOFF-1) and read_left_set==False:
            read_left_coor=current_coor_index-(current_read_index-READEND_CUTOFF+1)
            read_left_set=True
            
        ## Setup right read border
        if current_read_index>=(read_length-READEND_CUTOFF) and read_right_set==False:
            read_right_coor=current_coor_index-(current_read_index-(read_length-READEND_CUTOFF))
            read_right_set=True
        
        if current_coor_index>=search_coor and search_coor_found==False:
            if data_type=="N":
                SNV_base="N"
                filtered_out=True
                False_Positive_Type="splicing_skip"
                
            elif data_type=="D":
                SNV_base="D"
                #filtered_out=True
                False_Positive_Type="deletion"
                
            else:
                SNV_read_index=current_read_index-(current_coor_index-search_coor)
                SNV_base=read[SNV_read_index]
                SNV_quality=read_quality[SNV_read_index]
                if (ord(SNV_quality)-33) < 15:
                    filtered_out=True
                    False_Positive_Type="low_quality_call"
            search_coor_found=True
    
    ## check near splice junction   
    for border_coor in border_list:
        border_distance=abs(border_coor-search_coor)
        if border_distance<=(SPLICEJUNCTION_CUTOFF-1):
            False_Positive_Type="near_splice_junction"
            filter_out=True
            break
    
    ## check near read border
    if search_coor<=read_left_coor or search_coor>=read_right_coor:
        False_Positive_Type="near_read_border"
        filter_out=True
        
    if filter_out==False and SNV_base!="D":
        False_Positive_Type="NA"
        
    return filter_out,SNV_base,False_Positive_Type  

def retrieve_geneid_from_ANONNOVAR(tmp_gene):
    tmp_gene_list=[]
    
    ##First step resolve the potential quote issue 
    if tmp_gene.count('"')==2:
        gene=tmp_gene[1:-1]
    else:
        gene=tmp_gene
    
    ##2nd: resolve potential multiple gene issue 
    if gene.count(";")>0:
        gene_info_list=gene.split(";")
    else:
        gene_info_list=[gene]
    
    for gene_info in gene_info_list:
        if gene_info.count("(")==0 and gene_info.count(",")==0:
            ## this is case I: single gene id
            gene_id=gene_info
            if gene_id not in tmp_gene_list:
                tmp_gene_list.append(gene_id)
        elif gene_info.count("(")==0 and gene_info.count(",")>0:
            ## this is case II: multiple gene ids, separate by ,
            gene_id_list=gene_info.split(",")
            for gene_id in gene_id_list:
                if gene_id not in tmp_gene_list:
                    tmp_gene_list.append(gene_id)
        elif gene_info.count("(")>0 and gene_info.count(",")==0:
            ## this is case III: () additional info case
            gene_id_list=gene_info.split("(")
            gene_id=gene_id_list[0]
            if gene_id not in tmp_gene_list:
                tmp_gene_list.append(gene_id)
        elif gene_info.count("(")==1 and gene_info.count(",")>0:
            ## this is case IIIb
            gene_id_list=gene_info.split("(")
            gene_id=gene_id_list[0]
            if gene_id not in tmp_gene_list:
                tmp_gene_list.append(gene_id)
            gene_id_list=gene_info.split(")")
            if gene_id_list[1].count(",")>0:
                tmp_item=gene_id_list[1]
                tmp_item_list=tmp_item.split(",")
                gene_id=tmp_item_list[1]
                if gene_id not in tmp_gene_list:
                    tmp_gene_list.append(gene_id)
        else:
            ### this is case IIIc
            gene_info_list=gene_info.split("(")
            for item in gene_info_list:
                if item.count(")")==0:
                    gene_id=item
                    if gene_id not in tmp_gene_list:
                        tmp_gene_list.append(gene_id)
                else:
                    item_list=item.split(")")
                    new_item=item_list[1]
                    if new_item.count(",")==1:
                        gene_id=new_item[1:]
                        if gene_id not in tmp_gene_list:
                            tmp_gene_list.append(gene_id)
        
    return tmp_gene_list


def alt_freq(result,alt,ref):
    nt_dict=dict()
    for nt in nt_list:
        nt_dict[nt]=0
    
    result_list=result.split(":")
    total_count=0
    for nt_count_info in result_list:
        nt_count_info_list=nt_count_info.split("_")
        #print nt_count_info_list
        nt=nt_count_info_list[0]
        if nt!="D":
            count=int(nt_count_info_list[1])
            nt_dict[nt]=count
            total_count+=count
    
    alt_count=nt_dict[alt]
    if total_count!=0:
        alt_percentage=int(100.00*alt_count/total_count)
    else:
        alt_percentage=0
    
    return alt_percentage

def alt_freq_v2(result,alt,ref):
    nt_dict=dict()
    for nt in nt_list:
        nt_dict[nt]=0
    
    result_list=result.split(":")
    total_count=0
    for nt_count_info in result_list:
        nt_count_info_list=nt_count_info.split("_")
        #print nt_count_info_list
        nt=nt_count_info_list[0]
        if nt!="D":
            count=int(nt_count_info_list[1])
            nt_dict[nt]=count
            total_count+=count
    
    alt_count=nt_dict[alt]
    if total_count!=0:
        alt_percentage=int(100.00*alt_count/total_count)
    else:
        alt_percentage=-1
    
    return alt_percentage

def extract_coverage(data_info):
    coverage=0
    result_list=data_info.split(":")
    for nt_count_info in result_list:
        nt_count_info_list=nt_count_info.split("_")
        #print nt_count_info_list
        nt=nt_count_info_list[0]
        if nt!="D":
            count=int(nt_count_info_list[1])
            coverage+=count
    return coverage

def ExtractCount2Dict(data_info):
    count_dict = dict()
    result_list=data_info.split(":")
    for nt_count_info in result_list:
        nt_count_info_list=nt_count_info.split("_")
        #print nt_count_info_list
        nt=nt_count_info_list[0]
        if nt!="D":
            count=int(nt_count_info_list[1])
            count_dict[nt]=count
    return count_dict

def bam_index(bam_infile):
    '''
    function to test whether a index file exists for a bam file
    Status: working
    '''
    ## test if index exist or not
    bam_index_existence = False
    samtools_index = bam_infile+".bai"
    picard_index = bam_infile[:-1]+"i"
    if os.path.isfile(samtools_index) or os.path.isfile(picard_index):
        bam_index_existence = True
    return bam_index_existence

def fragment_length_SAMTOOLS(bam_infile):
    '''
    function: to obtain sequencing length from 
    Status: development
    '''
    
def extract_median_coverage(data_info_list):
    count_list=[]
    for data in data_info_list:
        sample_coverage=extract_coverage(data)
        count_list.append(sample_coverage)
    median_coverage=numpy.median(count_list)
    return median_coverage

def flag2fullbinary(flag):
    flag_binary = bin(int(flag))
    flag_binary_len = len(flag_binary)-2
    final_binary = "0"*(12-flag_binary_len)+flag_binary[2:]
    return final_binary

def flag_filter(flag_int,position):
    ## relative position
    ## for example: strand is position 5
    ## with 1 will be false, and filter out
    ## with 0, will be true and keep
    flag_binary = bin(flag_int)
    flag_length = len(flag_binary)-2
    if flag_length<position:
        return True
    else:
        flag_position=flag_binary[(-1)*position]
        if flag_position=="1":
            return False
        else:
            return True

def coverage_list2simple_dict(coverage_list,based_system=0):
    ##
    list_length=len(coverage_list)
    merged_start_list=[]
    merged_end_list=[]
    merged_dict=dict()
    for index in range(1,list_length+1):
        if index==1:
            previouse_coverage=coverage_list[index-1]
            index_start=index
            index_end= index
            pass
        else:
            current_coverage=coverage_list[index-1]
            if previouse_coverage==current_coverage:
                index_end=index
                pass
            else:
                if based_system==0:
                    final_index_start=index_start-1
                    final_index_end=index_end
                else:
                    final_index_start=index_start
                    final_index_end=index_end
                merged_start_list.append(final_index_start)
                merged_end_list.append(final_index_end)
                merged_dict[final_index_start]=previouse_coverage
                index_start=index
                index_end=index
                previouse_coverage=current_coverage
    return merged_start_list,merged_end_list,merged_dict


def coverage2bedGraph(coverage_dict,output_handle):
    for ref_ID in coverage_dict.keys():
        coverage_data_list = coverage_dict[ref_ID][1:]
        start_list,end_list,data_dict=coverage_list2simple_dict(coverage_data_list)
        for index in range(len(start_list)):
            start_index=start_list[index]
            coverage=data_dict[start_index]
            if coverage>0:
                output_handle.write(ref_ID+'\t')
                end_index=end_list[index]
                output_handle.write(str(start_index)+'\t')
                output_handle.write(str(end_index)+'\t')
                output_handle.write(str(coverage)+'\n')