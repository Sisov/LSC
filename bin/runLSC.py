#!/usr/bin/python

import argparse
import sys
import os
from numpy import * # Where is this being used?
import datetime
from re import * # Where is this being used?
from copy import * # Where is this being used?
import threading #This should be changed to a Pool later
import string
import commands # This is used for running some line counts later and can probably go away when those do.
from SequenceBasics import GenericFastaFileReader, GenericFastqFileReader
from multiprocessing import cpu_count
import subprocess
from random import randint

def log_print(print_str):
    os.system("echo '" + str(print_str) + "'")

def log_command(print_str):
    os.system("echo " + str(print_str))
    os.system(str(print_str))


################################################################################

def GetPathAndName(pathfilename):
    full_path = os.path.abspath(pathfilename)
    [path, filename] = os.path.split(full_path)
    path += '/'
    return path, filename

def Readcfgfile(cfg_filename):
    results = {}
    cfg = open(cfg_filename,'r')
    for line in cfg:
        line = line.strip()
        if line=='':
            continue
        if not line[0]=='#':
            ls = line.split('=')
            log_print(ls)
            if len(ls)>2:
                log_print('warning: too many = in cfg file')
            results[ls[0].strip()] = ls[1].strip()
    cfg.close()
    return results


def main():
  version = "2.alpha"
  parser = argparse.ArgumentParser(description="LSC "+version+": Correct errors (e.g. homopolymer errors) in long reads, using short read data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--long_reads',required=True,help="FASTAFILE Long reads to correct")
  parser.add_argument('--short_reads',nargs='+',help="FASTAFILE Short reads used to correct the long reads")
  parser.add_argument('--threads',type=int,default=0,help="Number of threads (Default = cpu_count)")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default='/tmp',help="FOLDERNAME where temporary files can be placed")
  group.add_argument('--specific_tempdir',help="FOLDERNAME of exactly where to place temproary folders")
  parser.add_argument('-o','--output',required=True,help="FOLDERNAME where output is to be written")
  parser.add_argument('--mode',default=0,type=int,help="0: run through")
  parser.add_argument('--sort_mem_max',type=int,help="-S option for memory in unix sort")
  args = parser.parse_args()
  if args.threads == 0:
    args.threads = cpu_count()

  run_pathfilename = os.path.realpath(__file__)

  # Establish running parameters here.  Some of these may need to become command line options.
  python_path = "/usr/bin/python"
  mode = args.mode
  LR_pathfilename = args.long_reads
  SR_pathfilename = args.short_reads[0]
  SR_filetype = 'fa'
  LR_filetype = 'fa'
  temp_foldername = 'temp'
  output_foldername = 'output'
  I_RemoveBothTails = "N"
  MinNumberofNonN = "40"
  MaxN = "1"
  I_nonredundant = "N"
  SCD = 20
  sort_max_mem = "-1"
  clean_up = 0
  max_error_rate = "12"
  bowtie2_options = "--end-to-end -a -f -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.08 --no-unal --omit-sec-seq"

  sys.stderr.write("=== Welcome to LSC " + version + " ===\n")
  
  ################################################################################
  # Make the temporary directory to work in
  # and output directory
  temp_foldername = args.tempdir.rstrip('/') + '/'
  if args.specific_tempdir:
    temp_foldername = args.specific_tempdir.rstrip('/')+'/'
  else:
    # We need to make a very temporary directory
    rname = "lsc."+str(randint(1,1000000000))
    if not os.path.isdir(temp_foldername+rname):
      os.makedirs(temp_foldername+rname)
      temp_foldername+=rname+'/'
  output_foldername = args.output.rstrip('/')+'/'

  if not os.path.isdir(output_foldername):
    os.makedirs(output_foldername)

  if not os.path.isdir(temp_foldername):
    if mode == 2:
      sys.stderr.write("Error: temp folder does not exist.\n")
      sys.stderr.write("Note: You need to run LSC in mode 1 or 0 before running in mode 2.\n")
      sys.exit()
    else: 
      os.makedirs(temp_foldername)
  if not os.path.isdir(temp_foldername + "log"):
    os.makedirs(temp_foldername+"log")

  bin_path, run_filename = GetPathAndName(run_pathfilename)
  LR_path, LR_filename = GetPathAndName(LR_pathfilename)
  SR_path, SR_filename = GetPathAndName(SR_pathfilename)

  python_bin_path = python_path + " " + bin_path

  t0 = datetime.datetime.now()

  # LSC does have an option for a 'clean_up' step with clean_up.py but I don't recall it being documented 
  # note that clean_up.py took as inputs the temp folder, and the two cpu_counts
        
  ################################################################################
  # Remove duplicate short reads first
    
  if mode == 0 or mode == 1:    
    if I_nonredundant == "N" and SR_filetype != "cps":  # If we go in, we want to get a unique set
        remove_duplicate_short_reads(SR_filetype,SR_pathfilename,temp_foldername,bin_path,args)
        sys.stderr.write(str(datetime.datetime.now()-t0))
        
  ##########################################
  # Split the SR FASTQ across CPUs
  ext_ls=[]
  for i in range(args.threads):
    ext_ls.append( '.' + string.lowercase[i / 26] + string.lowercase[i % 26] )

  if ((mode == 0) or 
    (mode == 1)):
    
        log_print("===split SR:===")    
        if (SR_filetype == "cps"):
            SR_NL = int(commands.getstatusoutput('wc -l ' + SR_pathfilename)[1].split()[0])
            Nsplitline = 1 + (SR_NL / args.threads)
            if (Nsplitline % 2 == 1):
                Nsplitline +=1
                    
            splitSR_cmd = "split -l " + str(Nsplitline) + " " + SR_pathfilename + " " + temp_foldername + "SR.fa."
            log_command(splitSR_cmd)
            for ext in ext_ls:
                mv_cmd = "mv " + temp_foldername + "SR.fa" + ext + " "  + temp_foldername + "SR.fa" + ext + ".cps"
                log_command(mv_cmd)
            splitSR_cmd = "split -l " + str(Nsplitline/2) + " " + SR_path + "SR.fa.idx " + temp_foldername + "SR.fa."
            log_command(splitSR_cmd)
            for ext in ext_ls:
                mv_cmd = "mv " + temp_foldername + "SR.fa" + ext + " "  + temp_foldername + "SR.fa" + ext + ".idx"
                log_command(mv_cmd)
                            
            log_print(str(datetime.datetime.now()-t0))
        else:
            SR_NL = int(commands.getstatusoutput('wc -l ' + SR_pathfilename)[1].split()[0])
            Nsplitline = 1 + (SR_NL / args.threads)
            if ( SR_filetype == "fa"):
                if (Nsplitline % 2 == 1):
                    Nsplitline +=1
            elif ( SR_filetype == "fq"):
                if (Nsplitline % 4 != 0):
                    Nsplitline += (4 - (Nsplitline % 4))
            else:
                log_print("Err: invalid filetype for short reads")
                exit(1)
            
            splitSR_cmd = "split -l " + str(Nsplitline) + " " + SR_pathfilename + " " + temp_foldername + "SR.fa."
            log_command(splitSR_cmd)
                        
            log_print(str(datetime.datetime.now()-t0))

  SR_filename = "SR.fa"
  if (SR_filetype == "cps"):
    SR_cps_pathfilename = SR_path +  "SR.fa"
  else:
    SR_cps_pathfilename = temp_foldername + "SR.fa"

  ##########################################
  # Run HC over the split SR input
  if ((mode == 0) or 
    (mode == 1)):

    if (SR_filetype != "cps"):
        log_print("===compress SR.??:===")    
        
        i = 0
        T_compress_SR_ls = []
        for ext in ext_ls:
            compress_SR_cmd = python_bin_path + "compress.py -MinNonN=" + MinNumberofNonN + " -MaxN=" + MaxN + " " + SR_filetype + " " + temp_foldername + SR_filename + ext + " " + temp_foldername + SR_filename + ext + "."
            T_compress_SR_ls.append( threading.Thread(target=log_command, args=(compress_SR_cmd,)) )
            T_compress_SR_ls[i].start()
            i += 1
        for T in T_compress_SR_ls:
            T.join()
        
        log_print(str(datetime.datetime.now()-t0))
        ####################
        # Remove temporary SR split files
        for ext in ext_ls:
            delSR_cmd = "rm " + temp_foldername + SR_filename + ext + " &"
            log_command(delSR_cmd)
        ####################

  ##########################################change output from compress.py and poolchr.py 
  # Remove the tails (shorter bits) from the LR, which SHOULD also be the overlapping DNA compliment
  if ((mode == 0) or 
    (mode == 1)):

    if I_RemoveBothTails == "Y":   
        log_print("===RemoveBothTails in LR:===")    
        RemoveBothTails_cmd = python_bin_path + "RemoveBothTails.py " + LR_filetype + " " + LR_pathfilename + " " + temp_foldername + "Notwotails_" + LR_filename 
        log_command(RemoveBothTails_cmd)
        LR_filetype = "fa"
        log_print(str(datetime.datetime.now()-t0))

    if I_RemoveBothTails == "Y":
        LR2fa_cmd = python_bin_path + "FASTA2fa.py " + temp_foldername + "Notwotails_" + LR_filename + " " + temp_foldername + "LR.fa"
        deltempLR_cmd = "rm " + temp_foldername + "Notwotails_" + LR_filename  
    else:
        if (LR_filetype == "fa"):
            LR2fa_cmd = python_bin_path + "FASTA2fa.py " + LR_pathfilename + " " + temp_foldername + "LR.fa"
        else:
            LR2fa_cmd = python_bin_path + "FASTQ2fa.py " + LR_pathfilename + " " + temp_foldername + "LR.fa"
            
    log_print(LR2fa_cmd)
    log_command(LR2fa_cmd)
    if I_RemoveBothTails == "Y":
        log_print(deltempLR_cmd)
        log_command(deltempLR_cmd)

  LR_filename = "LR.fa"

  ##########################################
  # Compress the long reads
  # Build the aligner index and then align short to long 
  if ((mode == 0) or 
    (mode == 1)):

    log_print(str(datetime.datetime.now()-t0))   
    
    log_print("===compress LR:===")    
    compress_LR_cmd = python_bin_path + "compress.py -MinNonN=" + MinNumberofNonN + " -MaxN=10000" + " fa " + temp_foldername + LR_filename + " " + temp_foldername + LR_filename +"."
    log_print(compress_LR_cmd)
    log_command(compress_LR_cmd)

    ####################
    LR_NL = int(commands.getstatusoutput('wc -l ' + temp_foldername + "LR.fa.cps")[1].split()[0])
    LR_NR = LR_NL / 2
    ####################
    delLR_cmd = "rm " + temp_foldername + "LR.fa"
    log_print(delLR_cmd)
    log_command(delLR_cmd)
    ####################

    log_print(str(datetime.datetime.now()-t0))

    # For now we can restrict ourselves to bowtie2.  hisat and bwa-mem are probably more likely to be used when this becomes modular
    log_print("===bowtie2 index LR:===")    
    bowtie2_index_cmd = "bowtie2-build -f " + temp_foldername + LR_filename + ".cps " + temp_foldername + LR_filename + ".cps"
    log_command(bowtie2_index_cmd)
        
    log_print(str(datetime.datetime.now()-t0))
        
        
    ##########################################
    log_print("===bowtie2 SR.??.cps:===")    
        
    i=0
    T_bowtie2_ls=[]
    for ext in ext_ls:
        bowtie2_cmd = "bowtie2 " + bowtie2_options + " -x " + temp_foldername + LR_filename + ".cps -U " + temp_foldername + SR_filename + ext + ".cps -S " + temp_foldername + SR_filename + ext + ".cps.sam" 
        T_bowtie2_ls.append( threading.Thread(target=log_command, args=(bowtie2_cmd,)) )
        T_bowtie2_ls[i].start()
        i+=1
    for T in T_bowtie2_ls:
        T.join()
    
    log_print(str(datetime.datetime.now()-t0))

    ##########################################
    # Convert the SAM file to a NAV file
    log_print("===samParser SR.??.cps.nav:===")
    
    i=0
    T_samParser_ls=[]
    for ext in ext_ls:
        samParser_cmd = (python_bin_path + "samParser.py " + temp_foldername + LR_filename + ".cps " + temp_foldername + SR_filename + ext + " " + 
                         temp_foldername + SR_filename + ext + ".cps.sam " + temp_foldername + SR_filename + ext + ".cps.nav " + max_error_rate + " ") 
        samParser_cmd += " T " 
        samParser_cmd += " > " + temp_foldername + SR_filename + ext + ".cps.samParser.log"
        T_samParser_ls.append( threading.Thread(target=log_command, args=(samParser_cmd,)) )
        T_samParser_ls[i].start()
        i+=1
    for T in T_samParser_ls:
        T.join()

    log_print(str(datetime.datetime.now()-t0))
    
  ##########################################
  # Build complete CPS and IDX files
    if (SR_filetype != "cps"):
        log_print("===cat SR.??.cps:===")    
        temp_filename_ls = []
        for ext in ext_ls:
            temp_filename_ls.append( temp_foldername + SR_filename + ext + ".cps" )
        log_command( "cat " + ' '.join(temp_filename_ls) + " > " + temp_foldername + SR_filename + ".cps" )
        log_print(str(datetime.datetime.now()-t0))
        ####################
        log_print("===cat SR.??.idx:===")    
        temp_filename_ls = []
        for ext in ext_ls:
            temp_filename_ls.append( temp_foldername + SR_filename + ext + ".idx" )
        log_command( "cat " + ' '.join(temp_filename_ls) + " > " + temp_foldername + SR_filename + ".idx" )
        log_print(str(datetime.datetime.now()-t0))
        ####################

    # Build complete SAM and NAV files
    ####################
    log_print("===cat SR.??.cps.sam :===")    
    temp_filename_ls = []
    for ext in ext_ls:
        temp_filename_ls.append( temp_foldername + SR_filename + ext + ".cps.sam" )
    log_command( "cat " + ' '.join(temp_filename_ls) + " > " + temp_foldername + SR_filename + ".cps.sam" )
    log_print(str(datetime.datetime.now()-t0))
    ####################
    log_print("===cat SR.??.cps.nav :===")    
    temp_filename_ls = []
    for ext in ext_ls:
        temp_filename_ls.append( temp_foldername + SR_filename + ext + ".cps.nav" )
    log_command( "cat " + ' '.join(temp_filename_ls) + " > " + temp_foldername + SR_filename + ".cps.nav" )
    log_print(str(datetime.datetime.now()-t0))
    ####################    

    ####################
    # Remove sam, alignment summary (.nav, .map) files per thread
    if (clean_up == 1):
        for ext in ext_ls:
            delSRnav_cmd = "rm " + temp_foldername + SR_filename + ext + ".cps.nav &"
            delSRsam_cmd = "rm " + temp_foldername + SR_filename + ext + ".cps.sam &"
            delSRidx_aa_cmd = "rm " + temp_foldername + SR_filename + ext + ".idx &"
            delSRcps_aa_cmd = "rm " + temp_foldername + SR_filename + ext + ".cps &"
            log_command(delSRcps_aa_cmd)
            log_command(delSRidx_aa_cmd)
            log_command(delSRnav_cmd)
            log_command(delSRsam_cmd)
            
    ####################
    for ext in ext_ls:
        log_command("mv " + temp_foldername + "SR.fa" + ext + ".cps.samParser.log " + temp_foldername + "log")
    ####################

  ####################

  # Return after alignment in case of mode 1
  if (mode == 1):
    exit(0)


  ##########################################
  log_print("===genLR_SRmapping SR.??.cps.nav:===")    
    
  genLR_SRmapping_cmd = python_bin_path + "genLR_SRmapping.py "  + temp_foldername + " " +  temp_foldername + SR_filename + ".cps.nav " + temp_foldername + LR_filename 
  genLR_SRmapping_cmd += " " + str(SCD) + " " + str(args.threads) + " " + str(sort_max_mem) 
  log_command(genLR_SRmapping_cmd)
  log_print(str(datetime.datetime.now()-t0))

  ##########################################
  log_print("===split LR_SR.map:===")    

  LR_SR_map_NR = int(commands.getstatusoutput('wc -l ' + temp_foldername +"LR_SR.map")[1].split()[0])

  if (LR_SR_map_NR == 0):
    log_print("Error: No short reads was aligned to long read. LSC could not correct any long read sequence.")
    exit(1)
    
  Nsplitline = 1 + (LR_SR_map_NR/args.threads)

  Nthread2_temp = int(LR_SR_map_NR)/int(Nsplitline)
  if ((LR_SR_map_NR % Nsplitline) != 0):
    Nthread2_temp += 1
  if (Nthread2_temp < args.threads):
    args.threads = Nthread2_temp

  ext2_ls=[]
  for i in range(args.threads):
    ext2_ls.append( '.' + string.lowercase[i / 26] + string.lowercase[i % 26] )
    
  splitLR_SR_map_cmd = "split -l " + str(Nsplitline) + " " + temp_foldername + "LR_SR.map" + ' ' + temp_foldername + "LR_SR.map" +"."
  log_command(splitLR_SR_map_cmd)

  log_print(str(datetime.datetime.now()-t0))

  ##########################################
  log_print("===correct.py LR_SR.map.??_tmp :===")    

  log_print(str(datetime.datetime.now()-t0))

  i=0
  T_correct_for_piece_ls=[]
  for ext in ext2_ls:
    correct_for_piece_cmd = python_bin_path + "correct_nonredundant.py " + temp_foldername + "LR_SR.map" + ext  + " " + temp_foldername + 'LR.fa.readname  > ' + temp_foldername + "LR_SR.map_emtry_ls" + ext
    T_correct_for_piece_ls.append( threading.Thread(target=log_command, args=(correct_for_piece_cmd,)) )
    T_correct_for_piece_ls[i].start()
    i+=1
  for T in T_correct_for_piece_ls:
    T.join()

  log_print(str(datetime.datetime.now()-t0))

  ####################
  if (clean_up == 1):
    for ext in ext2_ls:
        delLR_SR_map_aa_tmp_cmd = "rm " + temp_foldername + "LR_SR.map" + ext 
        log_command(delLR_SR_map_aa_tmp_cmd)
  ####################

  ##########################################

  log_print("===cat full_LR_SR.map.fa :===")    

  temp_filename_ls = []
  for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "full_LR_SR.map" + ext )
  log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "full_LR.fa" )

  log_print("===cat corrected_LR_SR.map.fa :===")    

  temp_filename_ls = []
  for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "corrected_LR_SR.map" + ext )
  log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "corrected_LR.fa" )

  log_print("===cat corrected_LR_SR.map.fq :===")    

  temp_filename_ls = []
  for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "corrected_LR_SR.map" + ext + '.fq' )
  log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "corrected_LR.fq" )

  log_print("===cat uncorrected_LR_SR.map.fa :===")    

  temp_filename_ls = []
  for ext in ext2_ls:
    temp_filename_ls.append( temp_foldername + "uncorrected_LR_SR.map" + ext  )
  log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldername + "uncorrected_LR.fa" )

  ####################
  if (clean_up == 1):
    for ext in ext2_ls:
        del_LR_SR_coverage_aa_cmd = "rm " + temp_foldername + "corrected_LR_SR.map" + ext + ".fq &"
        log_command(del_LR_SR_coverage_aa_cmd)
        delfull_LR_SR_map_aa_fa_cmd = "rm " + temp_foldername + "full_LR_SR.map" + ext 
        log_command(delfull_LR_SR_map_aa_fa_cmd)
        delcorr_LR_SR_map_aa_fa_cmd = "rm " + temp_foldername + "corrected_LR_SR.map" + ext
        log_command(delcorr_LR_SR_map_aa_fa_cmd)
        deluncorr_LR_SR_map_aa_fa_cmd = "rm " + temp_foldername + "uncorrected_LR_SR.map" + ext
        log_command(deluncorr_LR_SR_map_aa_fa_cmd)
  ####################

  ####################
  for ext in ext2_ls:
    log_command("mv " + temp_foldername + "LR_SR.map_emtry_ls" + ext + " " + temp_foldername + "log")
  ####################

# Pre: Requires the SR_filetype as fa or fq
#     SR_pathfilename as the loction of the SR file
#     temp_foldername is the location of the tempfolder
# Post: Writes SR_uniq.fa into tempfolder as a 

def remove_duplicate_short_reads(SR_filetype,SR_pathfilename,temp_foldername,bin_path,args):
  sys.stderr.write("=== sort and uniq SR data ===\n")
  if SR_filetype == 'fa':
    gfr = GenericFastaFileReader(SR_pathfilename)
  elif SR_filetype == 'fq':
    gfr = GenericFastqFileReader(SR_pathfilename)
  else:
    sys.stderr.write("ERROR: invalid SR_filetype "+SR_filetype+"\n")
    sys.exit()
  #Launch a pipe to store the unique fasta.  Its a little cumbersome but should
  #minimize memory usage.  Could go back and add the -S to sort if memory is a problem
  cmd = bin_path+"seq2uniqfasta.py --output "+temp_foldername+"SR_uniq.fa --tempdir "+temp_foldername
  if args.sort_mem_max:
    cmd += " --sort_mem "+str(args.sort_mem_max)
  p = subprocess.Popen(cmd.split(),stdin=subprocess.PIPE)
  while True:
    entry = gfr.read_entry()
    if not entry: break
    seq = entry['seq'].upper()
    p.stdin.write(seq+"\n")
  p.communicate()
  SR_pathfilename = temp_foldername + "SR_uniq.fa"
  return SR_pathfilename

if __name__=="__main__":
  main()
