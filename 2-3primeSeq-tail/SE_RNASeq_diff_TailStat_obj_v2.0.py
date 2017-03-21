# -*- coding: cp936 -*-
"""
The main purpose of this script is to generate stat for Single end Read
(Stranded) mainly CLIP-Seq reads regarding to transcript END.

=============================
Usage: python SE_RNASeq_diff_TailStat_obj_v0.1.alpha.py
-h help

-i files to processed                           *[No default value]

-b bed_infile to process                        *[No default value]

-R full reference length                        *[No default value]
     dict format 

-r GTF for transcription END file               *[No default value]
    (for differential END site)

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
1. *Tail.txt AND BED file Corresponding
   (Output Txt from SAM2SE_Tail_length_obj_v0.2.alpha.py) 

2. GTF file lable the transcription END for each reference (-r)
    differential END Site 
    format: first column reference gene/ID
    second column: the correct coordinate /range
    third column: the incorrect coordinate /range
    fourth column: the name of the gene/transcripts

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
6. No tail count 
<3> per gene stat file
<4> Per gene tail length stat

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
v0.0: 
"""

##Copyright
##By Liye Zhang
##Contact: bioliyezhang@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework

### Specific Functions definiation
def extract_info2coor_list(tail_info):
    region_count=tail_info.count(",")
    final_region_list = []
    if region_count>=1:
        region_list = tail_info.split(",")
        for region in region_list:
            if region.count("-")==1:
                start = int(region.split("-")[0])
                end = int(region.split("-")[1])
                for index in range(start,end+1):
                    final_region_list.append(index)
            else:
                final_region_list.append(int(region))
    else:
        if tail_info.count("-")==1:
            start = int(tail_info.split("-")[0])
            end = int(tail_info.split("-")[1])
            for index in range(start,end+1):
                final_region_list.append(index)
        else:
            final_region_list.append(int(tail_info))
    return final_region_list

def list2fastq(data_list,prefix,infile_obj):
    outfile_name=infile_obj.outputfilename_gen(prefix,"fq") ##create output file
    outfile_path=OUTPUT_PATH+"/"+outfile_name
    outfile_obj=GeneralFile_class(outfile_path)                           
    outfile_obj.output_handle_gen()
    for item in data_list:
        ID=item[0]
        seq=item[1]
        outfile_obj.handle.write("@"+ID+'\n')
        outfile_obj.handle.write(seq+"\n")
        outfile_obj.handle.write("+\n")
        outfile_obj.handle.write("P"*len(seq)+'\n')
    outfile_obj.handle.close()

def specific_function(infiles,bed_infiles,reference,full_reference):
    
    ##Section I: Generate the gene annotation dictionary
    
    cmd_records=record_command_line()  ##record the command line    
    
    ## Module 1: Process the full reference (dict format)
    infile_obj=GTF_File_class(full_reference)  
    infile_obj.SKIP_HEADER=1  ##unique ID length
    infile_reader=infile_obj.reader_gen()
    filtered_f_range_dict = dict()   ## record the region needs to care
    filtered_r_range_dict = dict() 
    ref_chro_list = []             ## reference list 
    length_dict = dict()           ## reference length dict
    for row in infile_reader:
        ref_ID = row[1][3:]
        full_length = int(row[2][3:])
        filtered_f_range_dict[ref_ID]=[0]*full_length ## f
        filtered_r_range_dict[ref_ID]=[0]*full_length ## r
        length_dict[ref_ID] = full_length
        

    ## Module 2: Process the gene Diff End information (not canonical format)
    ## Need to adjust this..!
    infile_obj=GeneralFile_class(reference)  ##create file obj(class)
    infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
    infile_reader=infile_obj.reader_gen()
    ref_gene_list = []
    for row in infile_reader:
        ref_ID = row[1]
        if ref_ID not in ref_chro_list:
            ref_chro_list.append(ref_ID)
        gene = row[0]
        if gene not in ref_gene_list:
            ref_gene_list.append(gene)
        strand= row[4]
        correct_end_coor_info = row[2]
        correct_end_coor_list = extract_info2coor_list(correct_end_coor_info)
        cryptic_end_coor_info = row[3]
        cryptic_end_coor_list = extract_info2coor_list(cryptic_end_coor_info)
        for index in correct_end_coor_list:
            if strand=="+":
                filtered_f_range_dict[ref_ID][index-1]=gene+"+"
            else:
                filtered_r_range_dict[ref_ID][index-1]=gene+"+"
        for index in cryptic_end_coor_list:
            if strand=="+":
                filtered_f_range_dict[ref_ID][index-1]=gene+"-"
            else:
                filtered_r_range_dict[ref_ID][index-1]=gene+"-"


    nt_list = ["A","T","G","C"]
    ## Module 3: Pair input file and BED file first
    if len(infiles)!=len(bed_infiles):
        print "Inconsistent Pairs,Please check output"
        sys.exit(0)
    else:
        paired_list = generate_paired_files_by_ID(infiles+bed_infiles,2,"txt","bed")
        ## return first file is txt , second file is bed

    ## Module 4: Data Processing
    for input_list in paired_list:
        infile = input_list[0]
        bed_infile = input_list[1]

        ## Module 4.1 : initiaze the dict that store the data
        print "Processing infile:", infile
        stat_dict_canonical = dict()
        stat_dict_cryptic = dict()
        len_dist_canonical_dict = dict()
        len_dist_cryptic_dict = dict() 
        tail_dict_canonical = dict()
        tail_dict_cryptic = dict()

        for gene in ref_gene_list:
            stat_dict_canonical[gene]=[0.0]*11 ## similar recording for each gene
            stat_dict_cryptic[gene]=[0.0]*11 ## similar recording for each gene
            len_dist_canonical_dict[gene] = []
            len_dist_cryptic_dict[gene] = []
            tail_dict_canonical[gene] = []
            tail_dict_cryptic[gene] = []
        
        #stat_dict[gene] = [0.0]*11 
        ## each reference: 
        ## within range count (1,index0)
        ## a Tail count (2,index1)
        ## u Tail count (3,index2) ; 
        ## Au count (4,index3) ; 
        ## unknown (5,index4) 
        ## nt (A,T,C,G) count 6-9 ; all atcg 10 
        ## no-tail count 11
        all_tail_count = 0 
        within_range_tail_count = 0 ## cryptic or canonical
        a_tail_count = 0 
        u_tail_count = 0
        au_tail_count = 0
        all_au_tail_count = 0
        tail_length_list = []
        nt_dict = dict()
        for nt in nt_list:
            nt_dict[nt]=0
        
        ## Module 4.2: include the notail bed file processing
        infile_obj=BEDFile_class(bed_infile)  ##create file obj(class)
        infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
        for row in infile_reader:
            chro=row[infile_obj.CHRO_COLUMN]
            coor=int(row[infile_obj.START_COLUMN])
            strand=row[infile_obj.STRAND_COLUMN]
            targeted_region=False
            if strand=="+" and filtered_f_range_dict[chro][coor]!=0:
                targeted_region=True
            elif strand=="-" and filtered_r_range_dict[chro][coor]!=0:
                targeted_region=True
            pass

            if targeted_region:
                if strand=="+":
                    tail_type=filtered_f_range_dict[chro][coor][-1]
                    gene = filtered_f_range_dict[chro][coor][:-1]
                else:
                    tail_type=filtered_r_range_dict[chro][coor][-1]
                    gene = filtered_r_range_dict[chro][coor][:-1]

                readcount = int(row[infile_obj.SCORE_COLUMN])
                within_range_tail_count+=readcount
                if tail_type=="+":
                    stat_dict_canonical[gene][10]+=readcount ## notail count
                    stat_dict_canonical[gene][0]+=readcount
                    len_dist_canonical_dict[gene]=len_dist_canonical_dict[gene]+[0]*readcount
                else:
                    stat_dict_cryptic[gene][10]+=readcount
                    stat_dict_cryptic[gene][0]+=readcount
                    len_dist_cryptic_dict[gene]=len_dist_cryptic_dict[gene]+[0]*readcount
                    
            ## to be continued:: test within the region add/update stat accordingly

        ## Module 4.3 : Process the Tail input file

        infile_obj=GeneralFile_class(infile)  ##create file obj(class)
        infile_obj.SKIP_HEADER=infile_skip    ##setup up the manual skip header if necessary
        infile_obj.SAMPLE_ID_LEN=unique_id_length  ##unique ID length
        infile_reader=infile_obj.reader_gen()  ##create the file reader to process infile
        for row in infile_reader:
            ref_ID = row[6]  
            if ref_ID.startswith("Tb"):
                if ref_ID.endswith("100nt")==False:
                    ref_ID=ref_ID+"_100nt" 
            all_tail_count+=1  
            strand = row[9]    
            #print "strand",strand 
            tail_seq = row[2]
            coor_end = int(row[8])
            tail_length = int(row[3])
            #final_coor = coor_end-tail_length
            if coor_end>length_dict[ref_ID]:
                coor_end= length_dict[ref_ID]

            targeted_region=False
            if strand=="+" and filtered_f_range_dict[ref_ID][coor_end-1]!=0:
                targeted_region=True
            elif strand=="-" and filtered_r_range_dict[ref_ID][coor_end-1]!=0:
                targeted_region=True
            else:
                pass

            if targeted_region:
                #print "come here"
                if strand=="+":
                    tail_type=filtered_f_range_dict[ref_ID][coor_end-1][-1]
                    gene = filtered_f_range_dict[ref_ID][coor_end-1][:-1]
                else:
                    tail_type=filtered_r_range_dict[ref_ID][coor_end-1][-1]
                    gene = filtered_r_range_dict[ref_ID][coor_end-1][:-1]

                within_range_tail_count+=1
                tail_ID = row[0]
                if tail_type=="+":
                    stat_dict_canonical[gene][0]+=1
                    len_dist_canonical_dict[gene].append(tail_length)
                    tail_dict_canonical[gene].append([tail_ID,tail_seq])
                else:
                    stat_dict_cryptic[gene][0]+=1
                    len_dist_cryptic_dict[gene].append(tail_length)
                    tail_dict_cryptic[gene].append([tail_ID,tail_seq])

                tail_length_list.append(tail_length)
                for nt in tail_seq:
                    nt_dict[nt]+=1

                    if tail_type=="+":
                        stat_dict_canonical[gene][9]+=1
                    else:
                        stat_dict_cryptic[gene][9]+=1

                    if nt=="A" or nt=="a":
                        if tail_type=="+":
                            stat_dict_canonical[gene][5]+=1
                        else:
                            stat_dict_cryptic[gene][5]+=1

                    elif nt=="T" or nt=="t":
                        if tail_type=="+":
                            stat_dict_canonical[gene][6]+=1
                        else:
                            stat_dict_cryptic[gene][6]+=1
                    
                    elif nt=="C" or nt=="c":
                        if tail_type=="+":
                            stat_dict_canonical[gene][7]+=1
                        else:
                            stat_dict_cryptic[gene][7]+=1
                    elif nt=="G" or nt=="g":
                        if tail_type=="+":
                            stat_dict_canonical[gene][8]+=1
                        else:
                            stat_dict_cryptic[gene][8]+=1
                    elif nt=="N" or nt=="n":
                        pass
                    else:
                        print "unrecognized nt",nt
                        sys.exit(0)
                
                ## check polyA tail
                tail_definition=row[-2]
                if tail_definition=="Atail":
                    if tail_type == "+":
                        stat_dict_canonical[gene][1]+=1
                    else:
                        stat_dict_cryptic[gene][1]+=1
                elif tail_definition=="Utail":
                    if tail_type == "+":
                        stat_dict_canonical[gene][2]+=1
                    else:
                        stat_dict_cryptic[gene][2]+=1
                elif tail_definition=="AUtail":
                    if tail_type == "+":
                        stat_dict_canonical[gene][3]+=1
                    else:
                        stat_dict_cryptic[gene][3]+=1
                elif tail_definition=="unknown":
                    if tail_type == "+":
                        stat_dict_canonical[gene][4]+=1
                    else:
                        stat_dict_cryptic[gene][4]+=1


                '''
                if tail_seq.count("A")>=int(float(tail_length)*0.90):
                    a_tail_count+=1
                    # a tail count (index1)
                    if tail_type=="+":
                        stat_dict_canonical[gene][1]+=1
                    else:
                        stat_dict_cryptic[gene][1]+=1
                if tail_seq.count("T")>=int(float(tail_length)*0.90):
                    u_tail_count+=1
                    # u tail count (index2)
                    if tail_type=="+":
                        stat_dict_canonical[gene][2]+=1
                    else:
                        stat_dict_cryptic[gene][2]+=1
                if (tail_seq.count("A")+tail_seq.count("T"))>=int(float(tail_length)*0.90):
                    au_tail_count+=1
                    # within range au tail count (index4)
                    if tail_type=="+":
                        stat_dict_canonical[gene][3]+=1
                    else:
                        stat_dict_cryptic[gene][3]+=1
                '''




            
        print "nt_dict",nt_dict
        ##Setup output file
        # outfile_name=infile_obj.outputfilename_gen("stats",OUTPUT_SUFFIX) ##create output file
        # outfile_path=OUTPUT_PATH+"/"+outfile_name
        # outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        # outfile_obj.RECORD=cmd_records                           
        # outfile_obj.output_handle_gen()    ##generate output handle       
        # outfile_obj.handle.write("#Item\tCount\n")
        # outfile_obj.handle.write("all tail count\t"+str(all_tail_count)+'\n')
        # outfile_obj.handle.write("within range tail count\t"+str(within_range_tail_count)+'\n')
        # within_range_tail_percent = float(within_range_tail_count)*100.00/float(all_tail_count)
        # outfile_obj.handle.write("within range tail percent\t"+str(within_range_tail_percent)+'\n')
        # outfile_obj.handle.write("A tail count\t"+str(a_tail_count)+'\n')
        # outfile_obj.handle.write("A tail percent within range tail\t")
        # a_tail_percent = 100.00 * float(a_tail_count)/float(within_range_tail_count)
        # outfile_obj.handle.write(str(a_tail_percent)+'\n')
        # outfile_obj.handle.write("U tail count\t"+str(u_tail_count)+'\n')
        # u_tail_percent = 100.00 * float(u_tail_count)/float(within_range_tail_count)
        # outfile_obj.handle.write(str(u_tail_percent)+'\n')
        # outfile_obj.handle.write("AU tail count\t"+str(au_tail_count)+'\n')
        # au_tail_percent = 100.00 * float(au_tail_count)/float(within_range_tail_count)
        # outfile_obj.handle.write(str(au_tail_percent)+'\n')
        # outfile_obj.handle.write("All + strand AU tail count\t"+str(all_au_tail_count)+'\n')
        #all_au_tail_percent = 100.00 * float(all_au_tail_count)/float(within_range_tail_count)
        #outfile_obj.handle.write(str(all_au_tail_percent)+'\n')
        # nt_all_count = float(sum(nt_dict.values()))
        #a_percent = 100.00 * float(nt_dict["A"])/nt_all_count
        #outfile_obj.handle.write("A %\t"+str(a_percent)+'\n')
        #c_percent = 100.00 * float(nt_dict["C"])/nt_all_count
        #outfile_obj.handle.write("C %\t"+str(c_percent)+'\n')
        #t_percent = 100.00 * float(nt_dict["T"])/nt_all_count
        #outfile_obj.handle.write("T %\t"+str(t_percent)+'\n')
        #g_percent = 100.00 * float(nt_dict["G"])/nt_all_count
        #outfile_obj.handle.write("G %\t"+str(g_percent)+'\n')
        # outfile_obj.handle.close()

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

        ## per gene output
        outfile_name=infile_obj.outputfilename_gen("per_gene_stat",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path) 
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()
        ## output header
        outfile_obj.handle.write("Description")
        for gene in ref_gene_list:
            ## output header
            outfile_obj.handle.write('\t'+gene+"_canonical")
            outfile_obj.handle.write('\t'+gene+"_cryptic")
        outfile_obj.handle.write('\n')

        description_list=["Within_range_tail_count","A tail count "
        ,"U tail count ","AU tail count","unknown tail count", "(within) A nt %",
        "U nt %","C nt %", "G nt %", "N nt","no_tail count"]

        for index in range(11):
            
            if index==0:
                outfile_obj.handle.write(description_list[index])
                for gene in ref_gene_list:
                    data=str(stat_dict_canonical[gene][index])
                    outfile_obj.handle.write('\t'+data)
                    data=str(stat_dict_cryptic[gene][index])
                    outfile_obj.handle.write('\t'+data)
                outfile_obj.handle.write('\n')
            if index in [1,2,3,4,10]:
                outfile_obj.handle.write(description_list[index])
                for gene in ref_gene_list:
                    data=str(stat_dict_canonical[gene][index])
                    outfile_obj.handle.write('\t'+data)
                    data=str(stat_dict_cryptic[gene][index])
                    outfile_obj.handle.write('\t'+data)
                outfile_obj.handle.write('\n')
                outfile_obj.handle.write(description_list[index]+"%")
                for gene in ref_gene_list:
                    if stat_dict_canonical[gene][0]==0:
                        total=1.00
                    else:
                        total=stat_dict_canonical[gene][0]
                    data=str(stat_dict_canonical[gene][index]*100.00/total)
                    outfile_obj.handle.write('\t'+data)
                    if stat_dict_cryptic[gene][0]==0:
                        total=1.00
                    else:
                        total=stat_dict_cryptic[gene][0]
                    data=str(stat_dict_cryptic[gene][index]*100.00/total)
                    outfile_obj.handle.write('\t'+data)
                outfile_obj.handle.write('\n')
            if index>=5 and index<=8:
                outfile_obj.handle.write(description_list[index])
                for gene in ref_gene_list:
                    if stat_dict_canonical[gene][9]==0:
                        total=1.00
                    else:
                        total=stat_dict_canonical[gene][9]
                    data=str(stat_dict_canonical[gene][index]*100.00/total)
                    outfile_obj.handle.write('\t'+data)
                    if stat_dict_cryptic[gene][9]==0:
                        total=1.00
                    else:
                        total=stat_dict_cryptic[gene][9]
                    data=str(stat_dict_cryptic[gene][index]*100.00/total)
                    outfile_obj.handle.write('\t'+data)
                outfile_obj.handle.write('\n')

        outfile_obj.handle.close()
        
        ## tail per gene length distribution
        outfile_name=infile_obj.outputfilename_gen("tail_length_dist_per_gene",OUTPUT_SUFFIX) ##create output file
        outfile_path=OUTPUT_PATH+"/"+outfile_name
        outfile_obj=GeneralFile_class(outfile_path)              ##create output obj
        outfile_obj.RECORD=cmd_records                           
        outfile_obj.output_handle_gen()    ##generate output handle
        outfile_obj.handle.write("#length")
        for gene in ref_gene_list:
            outfile_obj.handle.write('\t'+gene+"_canonical")
            outfile_obj.handle.write('\t'+gene+"_cryptic")
        outfile_obj.handle.write('\n')
        max_length= max(tail_length_list)
        for index in range(0,max_length+1):
            outfile_obj.handle.write(str(index))
            for gene in ref_gene_list:
                tail_length_list = len_dist_canonical_dict[gene]
                tail_count = tail_length_list.count(index)
                outfile_obj.handle.write('\t'+str(tail_count))
                tail_length_list = len_dist_cryptic_dict[gene]
                tail_count = tail_length_list.count(index)
                outfile_obj.handle.write('\t'+str(tail_count))
            outfile_obj.handle.write('\n')
        outfile_obj.handle.close()

        ## per gene output (fastq file)
        
        for gene in ref_gene_list:
            ## output cyptic
            prefix=gene+"-cyptic"
            if len(tail_dict_cryptic[gene])>0:
                list2fastq(tail_dict_cryptic[gene],prefix,infile_obj)
            prefix=gene+"-canonical"
            if len(tail_dict_canonical[gene])>0:
                list2fastq(tail_dict_canonical[gene],prefix,infile_obj)






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
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:b:s:S:r:R:u:d:D:j:I:t:p:L:o:O:z',["test="])
    for opt in optlist:
        if opt[0] == '-h':
            print __doc__; sys.exit(0)
        elif opt[0] == '-i': infile = opt[1]
        elif opt[0] == '-b': bed_infile = opt[1]
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
        bed_infiles=CurrentFolder_to_Infiles(INPUT_PATH, "bed")
    else:
        infiles=[infile]
        bed_infiles = [bed_infile]
    ##perform specific functions
    if reference!=None and full_reference!=None:
        specific_function(infiles,bed_infiles,reference,full_reference)
    
    
    
