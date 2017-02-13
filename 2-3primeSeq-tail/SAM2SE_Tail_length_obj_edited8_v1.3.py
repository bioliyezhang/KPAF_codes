# -*- coding: cp936 -*-
"""
The main purpose of this script is to extract Process Tail from single end Read.

=============================
Usage: python SAM2SE_Tail_length_obj_edited8_v1.2.py
-h help

-i files to processed                           *[No default value]

-r refrence dictionary file(.dict)              *[No default value]

-P primer (fasta format)                        *[No default value]

-R maximum Read length                          *[No default value]
   should be read length-1 

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
1.3:
fix incorrect handle of reverse strand 
1.2:
fix the H cigar string incorrect Tail retrieve and tail assignment
1.1:
fix over-estimate problem on ATail and underestimate Utail
due to incorrectly use round and int function

"""

##Copyright
##By Liye Zhang
##Contact: bioliyezhang@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework

### Specific Functions definiation
def detect_longest_shortU_tail(seq):

    U_count=0
    U_percentage_list=[]
    for index in range(len(seq)):
        nt=seq[index]
        if nt=="T":
            U_count+=1
        U_percentage = 100.00 * U_count / (index+1)
        if U_percentage<50.00:
            break
        else:
            U_percentage_list.append(U_percentage) 
    ##
    if len(U_percentage_list)==0:
        return 0
    else:
        max_U_percent = max(U_percentage_list)
        for index in range(len(U_percentage_list)):
            U_percentage=U_percentage_list[index]
            if U_percentage==max_U_percent:
                final_index=index+1
        return final_index

def tail_type_call(tail_seq):
    #tail_type="unknown"

    ## pattern definition
    #pattern_dict = dict()

    tail_length = len(tail_seq)
    U_count = tail_seq.count("T")
    if tail_seq.startswith("TT") and U_count>=round(0.90*tail_length):
        tail_type="Utail"
    elif tail_seq=="T":
        tail_type="Utail"
    else:
        ## (optional U + polyA) tail condition
        shortU_begin_len= detect_longest_shortU_tail(tail_seq)
        trim_Utail_tail_seq= tail_seq[shortU_begin_len:]
        trim_tail_len = len(trim_Utail_tail_seq)
        trim_A_count = trim_Utail_tail_seq.count("A")
        if trim_A_count>=round(0.90*trim_tail_len):
            tail_type="Atail"
        else:
            trim_U_count = trim_Utail_tail_seq.count("T")
            trim_AU_count = trim_A_count+trim_U_count
            if trim_AU_count>=round(0.80*trim_tail_len):
                tail_type="AUtail"
            else:
                tail_type="unknown"  
    
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


def specific_function(infiles,reference,primer_infile):
    ## infile: are SAM file
    ## reference: dict file generated for fasta file
    ## primer_infile: fasta for primer
    
    
    from numpy import median
    global pattern_list
    cmd_records=record_command_line()  ##record the command line 
    
    ## Part 1: Process the primer infile
    primer_dict = read_fasta(primer_infile)
    primer_list = primer_dict.values() 

    ## Part 2: Process the dict file for reference fasta
    infile_obj=GeneralFile_class(reference)  ##create file obj(class)
    infile_obj.SKIP_HEADER=1    ##setup up the manual skip header if necessary
    infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
    infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
    reference_length_dict= dict()
    for row in infile_reader:
        ID=row[1][3:]
        length=int(row[2][3:])
        reference_length_dict[ID]=length
    reference_list = reference_length_dict.keys()

    ## Part 3: Process SAM file
    for infile in infiles:

        print "Processing infile:", infile
        final_type_list= ["unknown","Utail","Atail","AUtail","notail",""]
        final_type_count_dict = dict()
        for tail_type in final_type_list:
            final_type_count_dict[tail_type]=0 ## number 

        ID_list = []
        tail_type_ID_dict = dict()
        ID_type_dict = dict() # Per ID record the Tail type
        for pattern in pattern_list:
            tail_type_ID_dict[pattern] = [] ##? not sure about this

        tail_data_dict = dict()
        notail_f_dict = dict()
        notail_r_dict = dict()
        for reference_ID in reference_list:
            length=reference_length_dict[reference_ID]
            notail_r_dict[reference_ID]=[0]*length ## record the no tail readcount
            notail_f_dict[reference_ID]=[0]*length

        
        ##Set up infile object
        infile_obj=SAM_File_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  
        
        ## round 1 Process (initilize the Dictionary)
        ## dictiontionary for tail length
        tail_length_dict =dict()
        full_seq_dict = dict() ## 
        for row in infile_reader:
            ID = row[infile_obj.QNAME_COLUMN]
            tail_length_dict[ID] = [1000,1000] # two value, one for maxcircle, one for edited
            ID_type_dict[ID] = "" ## make sure no empty value
            full_seq_dict[ID] = ""
        
        ## Round 2 record the different data (design to handle RPS12)
        infile_obj=SAM_File_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  
        for row in infile_reader:

            ID = row[infile_obj.QNAME_COLUMN]
            ref_ID = row[infile_obj.CHRO_COLUMN]
            
            ## design to capture which side
            if ref_ID=="max_circle":
                tail_index=0
            elif ref_ID.startswith("Tb"):
                tail_index=1
            else:
                tail_index=-1 ## non-maxcircle and rps12

            ## will make this a common function to handle flag
            flag = int(row[infile_obj.FLG_COLUMN])
            flag_binary = bin(flag)
            flag_binary_len = len(flag_binary)-2
            final_binary = "0"*(12-flag_binary_len)+flag_binary[2:]

            ## design on the potential tail length
            
            cigar_data = row[infile_obj.CIGAR_COLUMN]
            
            if cigar_data!="*":

                map_start_coor = int(row[infile_obj.COOR_COLUMN])
                cigar_extracted_list, read_length, match_length = read_cigar_coor_offset_withlist(cigar_data)
                min_length = len(row[infile_obj.SEQ_COLUMN])*MATCH_PERCENT_LIMIT/100

                ## only look into the Cigar that is not terribly complicated, current set to 10
                if len(cigar_extracted_list)<=CIGAR_UPPER_LIMIT and match_length>=min_length:
                    
                    full_seq = row[infile_obj.SEQ_COLUMN]
                    if final_binary[7]=="0":  ## 0 means forward, no reverse complement
                        direction="+"
                        if len(full_seq)>len(full_seq_dict[ID]):
                            full_seq_dict[ID]=full_seq
                        fragment_start = map_start_coor
                        fragment_end = map_start_coor + match_length -1
                        
                        first_cigar_string = cigar_extracted_list[-1]
                        if first_cigar_string[-1] in ["S","H"]:
                            tail_length = int(first_cigar_string[:-1])
                            if tail_index>=0:
                                if tail_length<tail_length_dict[ID][tail_index]:
                                    tail_length_dict[ID][tail_index]=tail_length
                                
                    else:
                        if len(full_seq)>len(full_seq_dict[ID]):
                            full_seq_dict[ID]=reverse_complementary(full_seq)
                        direction="-"
                        fragment_end = map_start_coor
                        fragment_start = map_start_coor+match_length-1
                        first_cigar_string = cigar_extracted_list[0]

                        if first_cigar_string[-1] in ["S","H"]:
                            if tail_index>=0:
                                tail_length = int(first_cigar_string[:-1])
                                if tail_length<tail_length_dict[ID][tail_index]:
                                    tail_length_dict[ID][tail_index]=tail_length 

        ## round 2 Process
        ##Set up infile object
        infile_obj=SAM_File_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  

        for row in infile_reader:
            ID = row[infile_obj.QNAME_COLUMN]
            ID_type_dict[ID] = "" ## what does this actually mean in here?
            ref_ID = row[infile_obj.CHRO_COLUMN]
            
            ## will make this a common function to handle flag
            flag = int(row[infile_obj.FLG_COLUMN])
            flag_binary = bin(flag)
            flag_binary_len = len(flag_binary)-2
            final_binary = "0"*(12-flag_binary_len)+flag_binary[2:]

            map_start_coor = int(row[infile_obj.COOR_COLUMN])
            cigar_data = row[infile_obj.CIGAR_COLUMN]
            if cigar_data!="*": 
                
                cigar_extracted_list, read_length, match_length = read_cigar_coor_offset_withlist(cigar_data)
                min_length = len(row[infile_obj.SEQ_COLUMN])*MATCH_PERCENT_LIMIT/100

                ## only look into the reads that are shorter than the total fragment
                if len(cigar_extracted_list)<=CIGAR_UPPER_LIMIT and match_length>=min_length:
                    start_with_designed_primer = False
                    
                    if len(ID_list)!=0:
                        if ID==ID_list[-1]:
                            pass
                        else:
                            ID_list.append(ID)         
                    else:
                        ID_list.append(ID)
                    
                    full_seq = full_seq_dict[ID]
                    #full_QC = row[infile_obj.QUAL_COLUMN]

                    tail_unmapped=False
                    if final_binary[7]=="0": 
                        ## forward strand
                        for primer_seq in primer_list:
                            if full_seq.startswith(primer_seq):
                                start_with_designed_primer=True
                                break
                        
                        direction="+"
                        fragment_start = map_start_coor
                        fragment_end = map_start_coor + match_length -1
                        
                        first_cigar_string = cigar_extracted_list[-1]
                        if first_cigar_string[-1] in ["S","H"]:
                            tail_length = int(first_cigar_string[:-1])
                            tail_unmapped=True
                            sense_tail_seq = full_seq[tail_length*(-1):]
                            #tail_QC = full_QC[tail_length*(-1):]
                            #QC_pass = QC_filter(sense_tail_seq,tail_QC)
                                
                    else:

                        direction="-"
                        fragment_end = map_start_coor
                        fragment_start = map_start_coor+match_length-1

                        first_cigar_string = cigar_extracted_list[0]
                        for primer_seq in primer_list:
                            if full_seq.startswith(primer_seq):
                                start_with_designed_primer=True
                                break


                        if first_cigar_string[-1] in ["S","H"]:
                            tail_length = int(first_cigar_string[:-1])
                            tail_unmapped=True
                            sense_tail_seq = full_seq[tail_length*(-1):]
                        

                    ## check the S,H count in cigar string
                    SH_fragment_count = count_SH_in_cigar(cigar_extracted_list)
                    #fragment_start = map_start_coor
                    if start_with_designed_primer==True:
                        if tail_unmapped==True and SH_fragment_count<3:
                            ## 
                            min_tail_length = min(tail_length_dict[ID])
                            if min_tail_length>0:
                                if len(sense_tail_seq)==min_tail_length:
                                    tail_type= tail_type_call(sense_tail_seq)
                                    ID_type_dict[ID] = tail_type
                                    tail_type_ID_dict[tail_type].append(ID)
                                    U_count = sense_tail_seq.count("T")
                                    A_count = sense_tail_seq.count("A")
                                    final_type_count_dict[tail_type]+=1
                                    tail_data_dict[ID] = [full_seq,sense_tail_seq,U_count,A_count,ref_ID,fragment_start,fragment_end,direction,tail_type]
                            else:
                                if read_length<=MAX_READLENGTH and first_cigar_string[-1]=="M":
                                    final_type_count_dict["notail"]+=1
                                    ID_type_dict[ID] = "notail"
                                    if direction=="+":
                                        notail_f_dict[ref_ID][fragment_end-1]+=1
                                    elif direction=="-":
                                        notail_r_dict[ref_ID][fragment_end-1]+=1
                                    else:
                                        print "direction unrecognized", direction
                                        sys.exit(0)
                        else:
                            ## with no tail read, require the READ length is smaller than the 300nt (max length)
                            
                            if read_length<=MAX_READLENGTH and first_cigar_string[-1]=="M":
                                final_type_count_dict["notail"]+=1
                                ID_type_dict[ID] = "notail"
                                if direction=="+":
                                    notail_f_dict[ref_ID][fragment_end-1]+=1
                                elif direction=="-":
                                    notail_r_dict[ref_ID][fragment_end-1]+=1
                                else:
                                    print "direction unrecognized", direction
                                    sys.exit(0)
                    else:
                        pass

        outfile_name=infile_obj.outputfilename_gen(prefix,OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle       
        outfile_obj.handle.write("#IDs\n")
        ID_set = set(ID_list)
        for seq_ID in ID_set:
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

        
        
        outfile_name=infile_obj.outputfilename_gen("Notail_Stat_BED6","bed") ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle  
        outfile_obj.handle.write("#ref_ID\tStart\tEND\tName\tReadCount\tStrand\n")
        ## output forward strand first
        for reference_ID in reference_list:
            length=reference_length_dict[reference_ID]
            for index in range(length):
                readcount=notail_f_dict[reference_ID][index]
                ID=reference_ID+"_"+str(index+1)+"_"+"p"
                outfile_obj.handle.write(reference_ID+'\t'+str(index)+'\t')
                outfile_obj.handle.write(str(index+1)+'\t'+ID+"\t"+str(readcount))
                outfile_obj.handle.write('\t+\n')
        ## output reverse strand first
        for reference_ID in reference_list:
            length=reference_length_dict[reference_ID]
            for index in range(length):
                readcount=notail_r_dict[reference_ID][index]
                ID=reference_ID+"_"+str(index+1)+"_"+"m"
                outfile_obj.handle.write(reference_ID+'\t'+str(index)+'\t')
                outfile_obj.handle.write(str(index+1)+'\t'+ID+"\t"+str(readcount))
                outfile_obj.handle.write('\t-\n')
        outfile_obj.handle.close()

        ##  generate the separate SAM file


        outfile_obj_dict = dict()
        for pattern in final_type_list:
            if pattern!="":
                outfile_name = infile_obj.outputfilename_gen(pattern,"sam")
            else:
                outfile_name = infile_obj.outputfilename_gen("unprocessed","sam")
            
            outfile_path = OUTPUT_PATH+"/"+outfile_name
            outfile_obj_dict[pattern] = GeneralFile_class(outfile_path) 
            outfile_obj_dict[pattern].output_handle_gen() 

        infile_obj=SAM_File_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_obj.AUTOSKIP_HEADER = False
        infile_reader=infile_obj.reader_gen()
        for row in infile_reader:
            if row[0][0]=="@":
                for pattern in final_type_list:
                    output_row(outfile_obj_dict[pattern].handle,row)
            else:
                ID = row[infile_obj.QNAME_COLUMN]
                ID_type = ID_type_dict[ID]
                output_row(outfile_obj_dict[ID_type].handle,row)

        for pattern in final_type_list:
            outfile_obj_dict[pattern].handle.close()

        ## 
        outfile_name=infile_obj.outputfilename_gen("stat",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle       
        outfile_obj.handle.write("#type\tCount\tPercentage\n")
        total_count=0
        for tail_type in final_type_list[:-1]:
            total_count+=final_type_count_dict[tail_type]

        for tail_type in final_type_list[:-1]:
            tail_count = final_type_count_dict[tail_type]
            tail_percent = 100.00 * tail_count/total_count
            outfile_obj.handle.write(tail_type+'\t'+str(tail_count)+'\t')
            outfile_obj.handle.write(str(tail_percent)+'\n')
        outfile_obj.handle.write("all\t"+str(total_count)+'\t100.00\n')
        outfile_obj.handle.close()








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
    primer_infile=None
    INPUT_PATH=os.getcwd()
    OUTPUT_PATH=os.getcwd()
    prefix="proper_PE_IDs"
    OUTPUT_SUFFIX="txt"
    CIGAR_UPPER_LIMIT=4
    MATCH_PERCENT_LIMIT=0.00
    MIN_SIZE = 0
    MAX_SIZE = 1000
    MAX_READLENGTH = 299
    reference = None

    global pattern_list
    global original_list
    # original_list = ["polyAG_30_pattern","polyAU_30_pattern","polyA_30_pattern",
    # "oligoAG_5_29_pattern","oligoAU_5_29_pattern","oligoA_5_29_pattern",
    # "polyU_5_200_pattern","shortA_2_4_pattern","shortU_2_4_pattern","polyUA_pattern"]
    # pattern_list = []
    # for pattern in original_list:
    #     pattern_list.append(pattern)
    #     ligation_type = pattern + "-ligation"
    #     pattern_list.append(ligation_type)
    pattern_list=["unknown","Utail","Atail","AUtail","notail",""]


    ###get arguments(parameters)
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:s:S:r:c:m:M:d:P:R:D:j:I:t:p:L:o:O:z',["test="])
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
        elif opt[0] == '-R': MAX_READLENGTH = int(opt[1])
        elif opt[0] == '-P': primer_infile = opt[1]
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

    
    specific_function(infiles,reference,primer_infile)

    
    
