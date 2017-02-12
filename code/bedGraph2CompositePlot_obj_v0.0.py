# -*- coding: cp936 -*-
"""
The main purpose of this script is to average the data from Bedgraph based on 
certain feature such as gene start and end, default handle -100 as the transcript
end site. 

=============================
Usage: python bedGraph2CompositePlot_obj_v0.0.py
-h help

-i files to processed                           *[No default value]

-r reference dict file (contain reference len info) *[No default value]

-f feature BED Format file                      [default:None] 

-o prefix of output files                       [default:"output"]

-s suffix for the files to be processed         [default value "bedGraph"]

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
1. BedGraphs (required)
2. Reference dict file (required)
3. Optional(feature BED Point format)

======================
output files:
1. data table (unsmoothed version)

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
    
def specific_function(infiles,reference,feature_infile):
    
    ##Section I: Generate the gene annotation dictionary
    
    cmd_records=record_command_line()  ##record the command line  

    ## Process Dictionary infile
    ref_len_dict = dict()
    ref_list = list()
    infile_obj=GeneralFile_class(reference)  ##create file obj(class)
    infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
    infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
    infile_reader=infile_obj.reader_gen()  
    for row in infile_reader:
        identifier=row[0]
        if identifier=="@SQ":
            ref_ID = row[1][3:]
            ref_length = int(row[2][3:])
            ref_len_dict[ref_ID] = ref_length
            ref_list.append(ref_ID)
    
    ## Process feature_infile
    feature_coor_dict=dict()
    if feature_infile==None:
        for ref_ID in ref_list:
            ref_length = ref_len_dict[ref_ID]
            feature_coor_dict[ref_ID] = (ref_length - 100)
    else:
        infile_obj=BEDFile_class(feature_infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  
        for row in infile_reader:
            ref_ID = row[infile_obj.CHRO_COLUMN] ## chromosome
            coor = int(row[infile_obj.END_COLUMN]) ## point coordinate the end one
            feature_coor_dict[ref_ID] = coor

    ## obtain the range for each side
    max_left = 0
    max_right = 0
    for ref_ID in ref_list:
        ref_length = ref_len_dict[ref_ID]
        ref_coor = feature_coor_dict[ref_ID]
        ref_left_coor = 1 - ref_coor
        ref_right_coor = ref_length - ref_coor
        if ref_left_coor<max_left:
            max_left = ref_left_coor
        if ref_right_coor>max_right:
            max_right = ref_right_coor

    ## Gene Count Per Coordinates
    gene_count_per_coor_dict = dict()
    coor_list = range(max_left,max_right+1)
    for coor in coor_list:
        gene_count_per_coor_dict[coor]=0
    for ref_ID in ref_list:
        ref_length = ref_len_dict[ref_ID]
        ref_coor = feature_coor_dict[ref_ID]
        ref_left_coor = 1 - ref_coor
        ref_right_coor = ref_length - ref_coor
        for coor in range(ref_left_coor,ref_right_coor+1):
            gene_count_per_coor_dict[coor]+=1

    ## Process the BedGraph input 
    ## Data structure: 
    ## Dictionary: the total read count for each coordinates
    ## Dictionary: the sample for each coordinates
    outfile_name="coor_geneCount.txt" ##create output file
    outfile_path=OUTPUT_PATH+"/"+outfile_name
    outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
    outfile_obj.RECORD=cmd_records                           
    outfile_obj.output_handle_gen()    ##generate output handle    
    for coor in coor_list:
        outfile_obj.handle.write(str(coor)+'\t')   
        outfile_obj.handle.write(str(gene_count_per_coor_dict[coor])+'\n')
    outfile_obj.handle.close()
    
    for infile in infiles:

        ## initlize read count
        readcount_per_coor_dict=dict()
        for coor in coor_list:
            readcount_per_coor_dict[coor]=0

        print "Processing infile:", infile
        ##Set up infile object
        infile_obj=BEDFile_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
        for row in infile_reader:
            ref_ID=row[infile_obj.CHRO_COLUMN]
            ref_coor = feature_coor_dict[ref_ID]
            start =int(row[infile_obj.START_COLUMN])
            end = int(row[infile_obj.END_COLUMN])
            start_index = start-ref_coor+1
            end_index = end-ref_coor
            if row[3].count("e")==0:
                readcount = int(row[3]) ## Specific for BedGraph format
            else:
                readcount = float(row[3])
            for coor_index in range(start_index,end_index+1):
                readcount_per_coor_dict[coor_index]+=readcount

            
        
        ##Setup output file
        outfile_name=infile_obj.outputfilename_gen("total_readcount",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle      
        for coor in coor_list:
            outfile_obj.handle.write(str(coor)+'\t')   
            outfile_obj.handle.write(str(readcount_per_coor_dict[coor])+'\n') 
        outfile_obj.handle.close()

        ##Setup output file
        outfile_name=infile_obj.outputfilename_gen("composite_average",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle  
        for coor in coor_list:
            outfile_obj.handle.write(str(coor)+'\t')   
            coor_average = float(readcount_per_coor_dict[coor])/float(gene_count_per_coor_dict[coor])
            outfile_obj.handle.write(str(coor_average)+'\n')      
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
    suffix="bedGraph"
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
    feature_infile=None
    
    ###get arguments(parameters)
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:s:S:f:r:d:D:j:I:t:p:L:o:O:z',["test="])
    for opt in optlist:
        if opt[0] == '-h':
            print __doc__; sys.exit(0)
        elif opt[0] == '-i': infile = opt[1]
        elif opt[0] == '-I': INPUT_PATH = opt[1]
        elif opt[0] == '-O': OUTPUT_PATH = opt[1]
        elif opt[0] == '-S': OUTPUT_SUFFIX = opt[1]
        elif opt[0] == '-s': suffix = opt[1]
        elif opt[0] == '-f': feature_infile = opt[1]
        elif opt[0] == '-d': sep_char =opt[1]
        elif opt[0] == '-D': sep_gene =opt[1]
        elif opt[0] == '-j': infile_skip= int(opt[1])
        elif opt[0] == '-r': reference = opt[1]
        elif opt[0] == '-o': prefix = opt[1]
        elif opt[0] == '-L': unique_id_length = int(opt[1])
        elif opt[0] == '--test': long_input = opt[1]
    
    #print "Test long input", long_input
    if infile==None:
        infiles=CurrentFolder_to_Infiles(INPUT_PATH, suffix)
    else:
        infiles=[infile]

    ##perform specific functions
    if reference!=None:
        specific_function(infiles,reference,feature_infile)
    else:
        print "reference dictionary file not provided"
        sys.exit(0)
    
    
    
