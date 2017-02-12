# -*- coding: cp936 -*-
"""
The main purpose of this script is to extract Process Tail from single end Read.

=============================
Usage: python SAM2SE_length_obj_v0.1.beta.py
-h help

-i files to processed                           *[No default value]

-r refrence dict file                         *[No default value]

-p minimum mapped percent of the read            [default: 0.00]

-m Minimum size range based on Bioanalyzer      [default: 0]

-M maximum size range based on Bioanalyzer      [default: 10000]

-c cigar string count                           [default:2]
    this means it allows cigar string to have equal or less that 2 fragments

-o prefix of output files                       [default:"precursor_length"]

-s suffix for the files to be processed         [default value "sam"]

-S output SUFFIX                                [default: "txt"]

-d sep_char within among columns                [default value '\t']

-j skip header lines                            [default value 1]
    if skip header lines = 0 indicates that there is no header line
    if skip header lines = n indicates that first n lines are header lines
    
-I input file path              [default: current folder]

-O output file path             [default: current folder]

-L unique_id length for the infile      [default value 2]

===================
input description:
input files:
1. SAM file

======================
output files:
1. overall precursor length
2. each transcript length , count

============================

Python & Module requirement:
Versions 2.x : 2.4 or above 
Module: No additional Python Module is required.

============================
Library file requirement:

Not Standalone version, few library file is required.

============================

External Tools requirement:

============================
command line example:

============================
versions update
0.1 beta:
fixed end coordinate error
"""

##Copyright
##By Liye Zhang
##Contact: bioliyezhang@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework

### Specific Functions definiation
def tail_type_call(tail_seq):
    tail_type="unknown"

    ## pattern definition
    pattern_dict = dict()

    pattern_dict["polyA_30_pattern"] = "^A{30,200}C{0,1}"
    pattern_dict["polyAG_30_pattern"] = "^A{30,200}G{1,200}"
    pattern_dict["polyAU_30_pattern"] = "^A{30,200}T{1,200}" 
    pattern_dict["oligoA_5_29_pattern"] = "^A{5,29}C{0,1}"
    pattern_dict["oligoAG_5_29_pattern"] = "^A{5,29}G{1,200}"
    pattern_dict["oligoAU_5_29_pattern"] = "^A{5,29}T{1,200}" 
    pattern_dict["polyU_5_200_pattern"] = "^T{5,200}"
    pattern_dict["shortA_2_4_pattern"] = "^A{2,4}$"
    pattern_dict["shortU_2_4_pattern"] = "^T{2,4}$"
    pattern_dict["polyUA_pattern"] = "^T{2,200}A{2,200}$"

    global original_list


    for pattern_ID in original_list:
        pattern = pattern_dict[pattern_ID]
        if re.search(pattern,tail_seq)!=None:
            left_seq = re.sub(pattern,"",tail_seq)
            if len(left_seq)==0:
                tail_type=pattern_ID
            else:
                tail_type=pattern_ID+"-ligation"
            break    
    
    return tail_type

def count_SH_in_cigar(cigar_data_list):
    SH_count = 0
    for cigar_data in cigar_data_list:
        if cigar_data[-1] in ["S","H"]:
            SH_count+=1
    return SH_count

def QC_filter(seq,qc_score,MIN_MEAN_QC=20,MAX_LOW_QC_COUNT=3,LOW_QC=10,MAX_N=0):
    result = True
    seq_length = len(seq)
    
    #module I: max N
    N_count = seq.count("N")+seq.count("n")
    if N_count>MAX_N:
        result=False

    #module II: Mean QC score
    qc_score_total=0
    qc_score_low_count=0
    for qc in qc_score:
        qc_score_number = ord(qc)-32
        qc_score_total+=qc_score_number
        if qc_score_number<=LOW_QC:
            qc_score_low_count+=1
    qc_score_mean = qc_score_total/seq_length

    if qc_score_mean<MIN_MEAN_QC:
        result=False
    if qc_score_low_count>=MAX_LOW_QC_COUNT:
        result=False
    return result


def specific_function(infiles,reference):
    
    ##Section I: Generate the gene annotation dictionary
    from numpy import median

    global pattern_list
    
    cmd_records=record_command_line()  ##record the command line 
    '''
    infile_obj=GTF_File_class(reference)  ##create file obj(class)
    infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
    infile_reader=infile_obj.reader_gen()
    ref_length_dict = dict()
    ref_chro_list = []
    for row in infile_reader:
        ref_ID = row[infile_obj.CHRO_COLUMN]
        ref_chro_list.append(ref_ID)
        length = int(row[infile_obj.END_COLUMN])
        ref_length_dict[ref_ID]=length
    '''
    
    for infile in infiles:
        print "Processing infile:", infile

        per_gene_count_dict = dict()d
        '''
        for ref_ID in ref_chro_list:
            per_gene_count_dict[ref_ID]=[0]*3 ## first one total tail count
            ## second one: AU count (90% AU,less than 2 C)
            ## tail (within in 20nt) count
        '''

        ID_list = []
        tail_type_ID_dict = dict()
        ID_type_dict = dict()
        for pattern in (pattern_list+["","unknown"]):
            tail_type_ID_dict[pattern] = []

        tail_data_dict = dict()
        
        ##Set up infile object
        infile_obj=SAM_File_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
        for row in infile_reader:
            #SAM_fragment_length = int(row[infile_obj.TLEN_COLUMN])
            ID = row[infile_obj.QNAME_COLUMN]
            ID_type_dict[ID] = ""
            ref_ID = row[infile_obj.CHRO_COLUMN]
            
            ## will make this a common function to handle flag
            flag = int(row[infile_obj.FLG_COLUMN])
            flag_binary = bin(flag)
            flag_binary_len = len(flag_binary)-2
            final_binary = "0"*(11-flag_binary_len)+flag_binary[2:]

            map_start_coor = int(row[infile_obj.COOR_COLUMN])
            cigar_data = row[infile_obj.CIGAR_COLUMN]
            if cigar_data!="*" and flag<2048:
                
                cigar_extracted_list, read_length, match_length = read_cigar_coor_offset_withlist(cigar_data)
                min_length = len(row[infile_obj.SEQ_COLUMN])*MATCH_PERCENT_LIMIT/100
                if len(cigar_extracted_list)<=CIGAR_UPPER_LIMIT and match_length>=min_length:
                    
                    
                    
                    if len(ID_list)!=0:
                        if ID==ID_list[-1]:
                            pass
                        else:
                            ID_list.append(ID)         
                    else:
                        ID_list.append(ID)
                    
                    full_seq = row[infile_obj.SEQ_COLUMN]
                    full_QC = row[infile_obj.QUAL_COLUMN]
                    ## check the read length
                    if final_binary[6]=="0": 
                        ## 
                        ## reverse complementary

                        
                        direction="+"
                        fragment_start = map_start_coor
                        fragment_end = map_start_coor + match_length -1
                        ## this is tail
                        
                        

                        tail_unmapped=False
                        
                        first_cigar_string = cigar_extracted_list[-1]
                        if first_cigar_string[-1] == "S":
                            
                            #tail_seq = full_seq[:int(first_cigar_string[:-1])]
                            #tail_QC = full_QC[:int(first_cigar_string[:-1])]
                            #tail_unmapped=True
                            #sense_tail_seq = reverse_complementary(tail_seq)
                            ## get the qC score checked.
                            #QC_pass = QC_filter(sense_tail_seq,tail_QC)
                            tail_length = int(first_cigar_string[:-1])
                            tail_unmapped=True
                            sense_tail_seq = full_seq[tail_length*(-1):]
                            tail_QC = full_QC[tail_length*(-1):]
                            QC_pass = QC_filter(sense_tail_seq,tail_QC)
                                
                    else:

                        direction="-"
                        fragment_end = map_start_coor
                        fragment_start = map_start_coor-match_length+1

                        first_cigar_string = cigar_extracted_list[0]
                        if first_cigar_string[-1] == "S":

                            tail_seq = full_seq[:int(first_cigar_string[:-1])]
                            tail_QC = full_QC[:int(first_cigar_string[:-1])]
                            tail_unmapped=True
                            sense_tail_seq = reverse_complementary(tail_seq)
                            ## get the qC score checked.
                            QC_pass = QC_filter(sense_tail_seq,tail_QC)
                        
                        
                        

                    ## check the S,H count in cigar string
                    SH_fragment_count = count_SH_in_cigar(cigar_extracted_list)
                    #fragment_start = map_start_coor

                    if tail_unmapped==True and QC_pass==True and SH_fragment_count<2:
                        ##
                        tail_type= tail_type_call(sense_tail_seq)
                        ID_type_dict[ID] = tail_type
                        tail_type_ID_dict[tail_type].append(ID)
                        U_count = sense_tail_seq.count("T")
                        A_count = sense_tail_seq.count("A")
                        tail_data_dict[ID] = [full_seq,sense_tail_seq,U_count,A_count,ref_ID,fragment_start,fragment_end,direction,tail_type]



        outfile_name=infile_obj.outputfilename_gen(prefix,OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle       
        outfile_obj.handle.write("#IDs\n")
        for seq_ID in ID_list:
            outfile_obj.handle.write(seq_ID+'\n')
        outfile_obj.handle.close()

        outfile_name=infile_obj.outputfilename_gen("Tail",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle       
        outfile_obj.handle.write("#fastq_ID\tfull_Rseq\ttail_seq\tTail_length\tUcount\t")
        outfile_obj.handle.write("Acount\tChro\tFragment_start\tFragment_end\tdirection\t")
        outfile_obj.handle.write("Type\n")
        
        for seq_ID in tail_data_dict.keys():
            #tail_count+=1
            outfile_obj.handle.write(seq_ID+'\t')
            data_list = tail_data_dict[seq_ID]
            outfile_obj.handle.write(data_list[0]+'\t')
            outfile_obj.handle.write(data_list[1]+'\t')
            outfile_obj.handle.write(str(len(data_list[1]))+'\t')
            outfile_obj.handle.write(str(data_list[2])+'\t') #U count
            outfile_obj.handle.write(str(data_list[3])+'\t') #A count
            outfile_obj.handle.write(str(data_list[4])+'\t') #Chromosome
            outfile_obj.handle.write(str(data_list[5])+'\t') #fragment start
            outfile_obj.handle.write(str(data_list[6])+'\t') #fragment end
            outfile_obj.handle.write(str(data_list[7])+'\t') #direction of transcript

            outfile_obj.handle.write(str(data_list[8])+'\n') #tail type
        outfile_obj.handle.close()

        ##  generate the output object and the handle
        # final_type_list= []
        # for pattern in (pattern_list+["unknown"]):
        #     if len(tail_type_ID_dict[pattern])>0:
        #         final_type_list.append(pattern)
        # final_type_list.append("")

        # outfile_obj_dict = dict()
        # for pattern in final_type_list:
        #     if pattern!="":
        #         outfile_name = infile_obj.outputfilename_gen(pattern,"sam")
        #     else:
        #         outfile_name = infile_obj.outputfilename_gen("unprocessed","sam")
            
        #     outfile_path = OUTPUT_PATH+"/"+outfile_name
        #     outfile_obj_dict[pattern] = GeneralFile_class(outfile_path) 
        #     outfile_obj_dict[pattern].output_handle_gen() 

        ## disable the SAM output
        # infile_obj=SAM_File_class(infile)  ##create file obj(class)
        # infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        # infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        # infile_obj.AUTOSKIP_HEADER = False
        # infile_reader=infile_obj.reader_gen()
        # for row in infile_reader:
        #     if row[0][0]=="@":
        #         for pattern in final_type_list:
        #             output_row(outfile_obj_dict[pattern].handle,row)
        #     else:
        #         ID = row[infile_obj.QNAME_COLUMN]
        #         ID_type = ID_type_dict[ID]
        #         output_row(outfile_obj_dict[ID_type].handle,row)

        # for pattern in final_type_list:
        #     outfile_obj_dict[pattern].handle.close()

        ## generate the Stat file for each file
        '''
        outfile_name=infile_obj.outputfilename_gen("Stat",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle  
        outfile_obj.handle.write("#tail_count\t"+str(tail_count)+'\n')
        '''








if __name__ == "__main__":
    ###Python General Module Import 
    import sys, csv, getopt, re
    import os
    import math
    from itertools import ifilter
    
    ##Liye own common function,class loading
    from Constant_Library import *
    from General_Library import *
    from File_Class import *  ###
    from Sequencing_Library import *
    
    OUTPUT_SEP_CHAR='\t'
    
    
            
                 
    #exit if not enough arguments
    if len(sys.argv) < 3:
        print __doc__
        sys.exit(0)
    
    ###set default value
    suffix="sam"
    infile=None
    infile_skip=0
    sep_char='\t'
    sep_gene=','
    header_file=None
    unique_id_length=2
    parameter_file=None
    INPUT_PATH=os.getcwd()
    OUTPUT_PATH=os.getcwd()
    prefix="proper_PE_IDs"
    OUTPUT_SUFFIX="txt"
    CIGAR_UPPER_LIMIT=4
    MATCH_PERCENT_LIMIT=0.00
    MIN_SIZE = 0
    MAX_SIZE = 1000
    reference = None

    global pattern_list
    global original_list
    original_list = ["polyAG_30_pattern","polyAU_30_pattern","polyA_30_pattern",
    "oligoAG_5_29_pattern","oligoAU_5_29_pattern","oligoA_5_29_pattern",
    "polyU_5_200_pattern","shortA_2_4_pattern","shortU_2_4_pattern","polyUA_pattern"]
    pattern_list = []
    for pattern in original_list:
        pattern_list.append(pattern)
        ligation_type = pattern + "-ligation"
        pattern_list.append(ligation_type)


    ###get arguments(parameters)
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:s:S:r:c:m:M:d:D:j:I:t:p:L:o:O:z',["test="])
    for opt in optlist:
        if opt[0] == '-h':
            print __doc__; sys.exit(0)
        elif opt[0] == '-i': infile = opt[1]
        elif opt[0] == '-I': INPUT_PATH = opt[1]
        elif opt[0] == '-O': OUTPUT_PATH = opt[1]
        elif opt[0] == '-S': OUTPUT_SUFFIX = opt[1]
        elif opt[0] == '-s': suffix = opt[1]
        elif opt[0] == '-d': sep_char =opt[1]
        elif opt[0] == '-r': reference =opt[1]
        elif opt[0] == '-m': MIN_SIZE = int(opt[1])
        elif opt[0] == '-M': MAX_SIZE = int(opt[1])
        elif opt[0] == '-p': MATCH_PERCENT_LIMIT = float(opt[1])
        elif opt[0] == '-c': CIGAR_UPPER_LIMIT = int(opt[1])
        elif opt[0] == '-j': infile_skip= int(opt[1])
        elif opt[0] == '-o': prefix = opt[1]
        elif opt[0] == '-L': unique_id_length = int(opt[1])
        elif opt[0] == '--test': long_input = opt[1]
    
    #print "Test long input", long_input
    if infile==None:
        infiles=CurrentFolder_to_Infiles(INPUT_PATH, suffix)
    else:
        infiles=[infile]

    
    specific_function(infiles,reference)

    
    
