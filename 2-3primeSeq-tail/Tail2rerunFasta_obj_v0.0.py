# -*- coding: cp936 -*-
"""
The main purpose of this script is to generate fasta file for Blat to run.
rerun the unknown and AUtail type into No U alignment based on Blat.

=============================
Usage: python Tail2runFasta_obj_v0.0.py
-h help

-i files to processed                           *[No default value]

-p parameter file               *[No default value]

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
1. Tail output txt

======================
output files:


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
    
def specific_function(infiles):
    
    ##Section I: Generate the gene annotation dictionary
    
    cmd_records=record_command_line()  ##record the command line    
        
    for infile in infiles:
        print "Processing infile:", infile
        ##Set up infile object
        infile_obj=GeneralFile_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
        
            
        
        ##Setup output file
        outfile_name=infile_obj.outputfilename_gen(prefix,OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        #outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle  
        for row in infile_reader:
            tail_type=row[-1]
            if tail_type in ["AUtail","unknown"]:
                ID=row[0]
                seq=row[1]
                outfile_obj.handle.write(">"+ID+'\n')
                for nt in seq:
                    if nt not in ["T","t"]:
                        outfile_obj.handle.write(nt)
                outfile_obj.handle.write("\n")     
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
    prefix="rerun"
    OUTPUT_SUFFIX="fa"
    
    ###get arguments(parameters)
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:s:S:r:d:D:j:I:t:p:L:o:O:z',["test="])
    for opt in optlist:
        if opt[0] == '-h':
            print __doc__; sys.exit(0)
        elif opt[0] == '-i': infile = opt[1]
        elif opt[0] == '-I': INPUT_PATH = opt[1]
        elif opt[0] == '-O': OUTPUT_PATH = opt[1]
        elif opt[0] == '-S': OUTPUT_SUFFIX = opt[1]
        elif opt[0] == '-s': suffix = opt[1]
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
    specific_function(infiles)
    
    
    
