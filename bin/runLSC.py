#!/usr/bin/python

import argparse
import sys
import os
import json
#from numpy import * # Where is this being used?
import datetime
from re import * # Where is this being used?
from copy import * # Where is this being used?
import threading #This should be changed to a Pool later
import string
import commands # This is used for running some line counts later and can probably go away when those do.
from SequenceBasics import GenericFastaFileReader, GenericFastqFileReader
from multiprocessing import cpu_count, Pool
import subprocess
from random import randint
from SequenceCompressionBasics import HomopolymerCompressionFactory
from SamToNavBasics import SamToNavFactory

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
  parser.add_argument('--minNumberofNonN',type=int,default=40,help="Minimum number of non-N characters in the compressed read")
  parser.add_argument('--maxN',type=int,help="Maximum number of Ns in the compressed read")
  parser.add_argument('--error_rate_threshold',type=int,default=12,help="Maximum percent of errors in a read to use the alignment")
  parser.add_argument('--short_read_coverage_threshold',type=int,default=20,help="Minimum short read coverage to do correction")
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
  I_nonredundant = "N"
  sort_max_mem = -1
  clean_up = 0
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
        sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")
        
  ############################################################################
  # Run compression over short reads
  if mode == 0 or mode == 1 and SR_filetype != "cps":
    sys.stderr.write("===compress SR:===\n")    
    # At this point we should have uniq fasta formatted reads
    # We may want to consider supportin a unique set of fastq reads also
    global SR_cps_fh
    global SR_idx_fh
    SR_cps_fh = open(temp_foldername+'SR.fa.cps','w')
    SR_idx_fh = open(temp_foldername+'SR.fa.idx','w')
    compress_SR(temp_foldername+'SR_uniq.fa',args)
    SR_cps_fh.close()
    SR_idx_fh.close()
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")
        
  ##########################################
  # Change the names of the long read files and save the names
  # Previously at this step LSC had an option to remove fragments of long reads
  # That were usually the start and ending tails
  # That step is for pacbio reads and could be done without short reads, 
  # and thus seems to fall out of the scope of the purpose of LSC
  if mode == 0 or mode == 1:
    sys.stderr.write("===rename LR:===\n")    
    gfr = GenericFastaFileReader(args.long_reads)
    of_lr = open(temp_foldername+'LR.fa','w')
    of_lr_readnames = open(temp_foldername+'LR.fa.readnames','w')
    z = 0
    while True:
      entry = gfr.read_entry()
      if not entry: break
      z += 1
      of_lr.write(">"+entry['name']+"\n"+entry['seq'].upper()+"\n")
      of_lr_readnames.write(str(z)+"\t"+entry['name']+"\n")
    of_lr.close()
    of_lr_readnames.close()
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")

  LR_filename = "LR.fa"

  ##########################################
  # Compress the long reads
  # Build the aligner index and then align short to long 
  # This is the step where parallelization aught to come in to play
  # so I am ending mode 1 at this step.
  if mode == 0 or mode == 2:
    sys.stderr.write("===compress LR, samParser LR.??.cps.nav:===\n")    
    total_batches = execute_LR(temp_foldername+'LR.fa',args,temp_foldername)

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

def compress_SR(short_read_file,args):
  gfr = GenericFastaFileReader(short_read_file)
  batchsize = 100000
  batch = []
  cpus = args.threads
  p = Pool(processes=cpus)
  while True:
    entry = gfr.read_entry()
    if not entry: break
    batch.append([entry['name'],entry['seq']])
    if len(batch) >= batchsize: # we have enough to process
      p.apply_async(compress_batch,args=(json.dumps(batch),args.minNumberofNonN,args.maxN,),callback=collect_results)
      #results = compress_batch(batch)
      batch = []
  if len(batch) >= 0:
    p.apply_async(compress_batch,args=(json.dumps(batch),args.minNumberofNonN,args.maxN,),callback=collect_results)
    #results = compress_batch(batch)
  p.close()
  p.join()
  gfr.close()

def collect_results(results):
  global SR_cps_fh
  global SR_idx_fh    
  for entry in results:
    SR_cps_fh.write(">"+entry[0]+"\n"+entry[1]+"\n")
    SR_idx_fh.write(entry[0]+"\t"+entry[2]+"\t"+entry[3]+"\n")
  return

def compress_batch(seq_batch_json,minNumberofNonN,maxN):
  seq_batch = json.loads(seq_batch_json)
  result = []
  hpcf = HomopolymerCompressionFactory()
  if minNumberofNonN: hpcf.set_MinNonN(minNumberofNonN)
  if maxN: hpcf.set_MaxN(maxN)
  for entry in seq_batch:
    res = hpcf.compress(entry[1])
    if not res: continue
    (cps, pos_ls, len_ls) = res
    result.append([entry[0],cps,pos_ls,len_ls])
  return result

def execute_LR(lr_filename,args,temp_foldername):
  gfr = GenericFastaFileReader(lr_filename)
  batchsize = 10000
  batch = []
  batch_number = 0
  while True:
    entry = gfr.read_entry()
    if not entry: break
    batch.append([entry['name'],entry['seq']])
    if len(batch) >= batchsize:
      batch_number += 1
      execute_batch(json.dumps(batch),batch_number,args.minNumberofNonN,args.maxN,temp_foldername,args.threads,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max)
      batch = []
  if len(batch) > 0:
    batch_number += 1
    execute_batch(json.dumps(batch),batch_number,args.minNumberofNonN,args.maxN,temp_foldername,args.threads,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max)
  sys.stderr.write("\n")
  return batch_number

def execute_batch(batch_json,batch_number,minNumberofNonN,maxN,temp_foldername,threads,error_rate_threshold,short_read_coverage_threshold,sort_max_mem):
  batch = json.loads(batch_json)
  #step 1 compression
  sys.stderr.write("... executing batch "+str(batch_number)+".1\r")
  hpcf = HomopolymerCompressionFactory()
  if minNumberofNonN: hpcf.set_MinNonN(minNumberofNonN)
  if maxN: hpcf.set_MaxN(maxN)
  of_cps = open(temp_foldername+"LR.fa."+str(batch_number)+'.cps','w')
  of_idx = open(temp_foldername+"LR.fa."+str(batch_number)+'.idx','w')
  for entry in batch:
    res = hpcf.compress(entry[1])
    if not res: continue
    (cps, pos_ls, len_ls) = res
    of_cps.write(">"+entry[0]+"\n"+cps+"\n")
    of_idx.write(entry[0]+"\t"+pos_ls+"\t"+len_ls+"\n")
  of_cps.close()
  of_idx.close()
  #step 2 build an index
  sys.stderr.write("... executing batch "+str(batch_number)+".2\r")
  of_log = open(temp_foldername+"LR.fa."+str(batch_number)+'.cps'+'.log','w')
  cmd1 = 'hisat-build '+temp_foldername+'LR.fa.'+str(batch_number)+'.cps '+temp_foldername+'LR.fa.'+str(batch_number)+'.cps.hisat-index'
  subprocess.call(cmd1.split(),stderr=of_log,stdout=of_log)
  #step 3 map short reads
  sys.stderr.write("... executing batch "+str(batch_number)+".3\r")
  cmd2 = 'hisat --reorder -f -p '+str(threads)+' -x '+temp_foldername+'LR.fa.'+str(batch_number)+'.cps.hisat-index -U '+temp_foldername+'SR.fa.cps -S '+temp_foldername+'LR.fa.'+str(batch_number)+'.cps.sam'
  subprocess.call(cmd2.split(),stderr=of_log,stdout=of_log)
  #step 4 convert sam to nav
  sys.stderr.write("... executing batch "+str(batch_number)+".4\r")
  of = open(temp_foldername+'LR.fa.'+str(batch_number)+'.cps.nav','w')
  stnf = SamToNavFactory()
  stnf.set_error_rate_threshold(error_rate_threshold)
  stnf.initialize_compressed_query(temp_foldername+'SR.fa.cps',temp_foldername+'SR.fa.idx')
  stnf.initialize_target(temp_foldername+'LR.fa.'+str(batch_number)+'.cps')
  with open(temp_foldername+'LR.fa.'+str(batch_number)+'.cps.sam') as inf:
    for line in inf:
      navs = stnf.sam_to_nav(line)
      if not navs: continue
      for oline in navs:
        of.write(oline+"\n")
  stnf.close()
  of.close()

  #step 5 nav to mapping
  python_bin_path = "python ../LSC/bin/"
  genLR_SRmapping_cmd = python_bin_path + "genLR_SRmapping.py "  + temp_foldername + " " +  temp_foldername + "LR.fa."+str(batch_number)+".cps.nav " + temp_foldername + 'LR.fa.'+str(batch_number)
  genLR_SRmapping_cmd += " " + str(short_read_coverage_threshold) + " " + str(threads) + " " + str(sort_max_mem) +" "
  genLR_SRmapping_cmd += temp_foldername+"LR.fa.readnames "+temp_foldername+"LR_SR.map."+str(batch_number)
  print genLR_SRmapping_cmd  
  of_log.close()

if __name__=="__main__":
  main()
