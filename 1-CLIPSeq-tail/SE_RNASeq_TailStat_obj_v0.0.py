# -*- coding: cp936 -*-
"""
The main purpose of this script is to generate stat for Single end Read
(Stranded) mainly CLIP-Seq reads regarding to transcript END.

=============================
Usage: python SE_RNASeq_TailStat_obj_v0.0.py
-h help

-i files to processed                           *[No default value]

-R full reference length                        *[No default value]

-r GTF for transcription END file               *[No default value]

-u upstream nt to allow                         [default: -100]

-d downstream nt to allow                       [default: 100]

-o prefix of output files                       [default:"output"]

-s suffix for the files to be processed         [default value "txt"]

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
1. Tail (Output Txt from SAM2SE_Tail_length_obj_v0.1.alpha.py)
2. GTF file lable the transcription END for each reference (-r)
3. Dict file (-R)

======================
output files:
<1> stat file
<2> Tail length stat
1. all count, the distribution of tail length (total reads), total count
2. the distribution of A,T,C,G 4 nucleotide for all tails 
3. the read count of Atail (>90 percent A), percentage
4. the read count of Utail (>90 percent T), percentage
5. UA tail (>90 percent AU)

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

"""

##Copyright
##By Liye Zhang
##Contact: bioliyezhang@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework

### Specific Functions definiation
    
def specific_function(infiles,reference,full_reference):
    
    ##Section I: Generate the gene annotation dictionary
    
    cmd_records=record_command_line()  ##record the command line    
    ## Process the full reference
    infile_obj=GTF_File_class(full_reference)  ##create file obj(class)
    infile_obj.SKIP_HEADER=1  ##unique ID length
    infile_reader=infile_obj.reader_gen()
    filtered_range_dict = dict()
    ref_chro_list = []
    length_dict = dict()
    for row in infile_reader:
        ref_ID = row[1][3:]
        full_length = int(row[2][3:])
        filtered_range_dict[ref_ID]=[0]*full_length
        length_dict[ref_ID] = full_length

    ##
    infile_obj=GTF_File_class(reference)  ##create file obj(class)
    infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
    infile_reader=infile_obj.reader_gen()
    
    for row in infile_reader:
        ref_ID = row[infile_obj.CHRO_COLUMN]
        ref_chro_list.append(ref_ID)
        length = int(row[infile_obj.END_COLUMN])
        for index in range(length+UPSTREAM_BORDER-1,length+DOWNSTREAM_BORDER):
            filtered_range_dict[ref_ID][index]=1

    nt_list = ["A","T","G","C"]
    for infile in infiles:
        print "Processing infile:", infile
        all_tail_count = 0
        within_range_tail_count = 0
        a_tail_count = 0
        u_tail_count = 0
        au_tail_count = 0
        all_au_tail_count = 0
        tail_length_list = []
        nt_dict = dict()
        for nt in nt_list:
            nt_dict[nt]=0
        ##Set up infile object
        infile_obj=GeneralFile_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
        for row in infile_reader:
            all_tail_count+=1
            strand = row[-2]
            if strand == "+":
                ref_ID = row[6]
                tail_seq = row[2]
                coor_end = int(row[8])
                tail_length = int(row[3])
                final_coor = coor_end-tail_length
                if final_coor>length_dict[ref_ID]:
                    final_coor= length_dict[ref_ID]

                if (tail_seq.count("A")+tail_seq.count("T"))>=int(float(tail_length)*0.90):
                    all_au_tail_count+=1

                if filtered_range_dict[ref_ID][final_coor-1]==1:
                    within_range_tail_count+=1
                    tail_length_list.append(tail_length)
                    for nt in tail_seq:
                        nt_dict[nt]+=1
                    ## check polyA tail
                    if tail_seq.count("A")>=int(float(tail_length)*0.90):
                        a_tail_count+=1
                    if tail_seq.count("T")>=int(float(tail_length)*0.90):
                        u_tail_count+=1
                    if (tail_seq.count("A")+tail_seq.count("T"))>=int(float(tail_length)*0.90):
                        au_tail_count+=1





            
        
        ##Setup output file
        outfile_name=infile_obj.outputfilename_gen("stats",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle       
        outfile_obj.handle.write("#Item\tCount\n")
        outfile_obj.handle.write("all tail count\t"+str(all_tail_count)+'\n')
        outfile_obj.handle.write("within range tail count\t"+str(within_range_tail_count)+'\n')
        within_range_tail_percent = float(within_range_tail_count)*100.00/float(all_tail_count)
        outfile_obj.handle.write("within range tail percent\t"+str(within_range_tail_percent)+'\n')
        outfile_obj.handle.write("A tail count\t"+str(a_tail_count)+'\n')
        outfile_obj.handle.write("A tail percent within range tail\t")
        a_tail_percent = 100.00 * float(a_tail_count)/float(within_range_tail_count)
        outfile_obj.handle.write(str(a_tail_percent)+'\n')
        outfile_obj.handle.write("U tail count\t"+str(u_tail_count)+'\n')
        u_tail_percent = 100.00 * float(u_tail_count)/float(within_range_tail_count)
        outfile_obj.handle.write(str(u_tail_percent)+'\n')
        outfile_obj.handle.write("AU tail count\t"+str(au_tail_count)+'\n')
        au_tail_percent = 100.00 * float(au_tail_count)/float(within_range_tail_count)
        outfile_obj.handle.write(str(au_tail_percent)+'\n')
        outfile_obj.handle.write("All + strand AU tail count\t"+str(all_au_tail_count)+'\n')
        all_au_tail_percent = 100.00 * float(all_au_tail_count)/float(within_range_tail_count)
        outfile_obj.handle.write(str(all_au_tail_percent)+'\n')
        nt_all_count = float(sum(nt_dict.values()))
        a_percent = 100.00 * float(nt_dict["A"])/nt_all_count
        outfile_obj.handle.write("A %\t"+str(a_percent)+'\n')
        c_percent = 100.00 * float(nt_dict["C"])/nt_all_count
        outfile_obj.handle.write("C %\t"+str(c_percent)+'\n')
        t_percent = 100.00 * float(nt_dict["T"])/nt_all_count
        outfile_obj.handle.write("T %\t"+str(t_percent)+'\n')
        g_percent = 100.00 * float(nt_dict["G"])/nt_all_count
        outfile_obj.handle.write("G %\t"+str(g_percent)+'\n')
        outfile_obj.handle.close()

        ## tail length distribution
        outfile_name=infile_obj.outputfilename_gen("tail_length_dist",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle
        max_length= max(tail_length_list)
        for index in range(1,max_length+1):
            tail_count = tail_length_list.count(index)
            outfile_obj.handle.write(str(index)+'\t'+str(tail_count)+'\n')
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
    suffix="txt"
    infile=None
    infile_skip=0
    sep_char='\t'
    sep_gene=','
    header_file=None
    unique_id_length=2
    parameter_file=None
    INPUT_PATH=os.getcwd()
    OUTPUT_PATH=os.getcwd()
    prefix="output"
    OUTPUT_SUFFIX="txt"
    reference=None
    full_reference=None
    UPSTREAM_BORDER = -100
    DOWNSTREAM_BORDER = 100
    
    ###get arguments(parameters)
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:s:S:r:R:u:d:D:j:I:t:p:L:o:O:z',["test="])
    for opt in optlist:
        if opt[0] == '-h':
            print __doc__; sys.exit(0)
        elif opt[0] == '-i': infile = opt[1]
        elif opt[0] == '-I': INPUT_PATH = opt[1]
        elif opt[0] == '-O': OUTPUT_PATH = opt[1]
        elif opt[0] == '-S': OUTPUT_SUFFIX = opt[1]
        elif opt[0] == '-s': suffix = opt[1]
        elif opt[0] == '-d': DOWNSTREAM_BORDER = int(opt[1])
        elif opt[0] == '-u': UPSTREAM_BORDER = int(opt[1])
        elif opt[0] == '-D': sep_gene =opt[1]
        elif opt[0] == '-j': infile_skip= int(opt[1])
        elif opt[0] == '-r': reference = opt[1]
        elif opt[0] == '-R': full_reference = opt[1]
        elif opt[0] == '-o': prefix = opt[1]
        elif opt[0] == '-L': unique_id_length = int(opt[1])
        elif opt[0] == '--test': long_input = opt[1]
    
    #print "Test long input", long_input
    if infile==None:
        infiles=CurrentFolder_to_Infiles(INPUT_PATH, suffix)
    else:
        infiles=[infile]

    ##perform specific functions
    if reference!=None and full_reference!=None:
        specific_function(infiles,reference,full_reference)
    
    
    
