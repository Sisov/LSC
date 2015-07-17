#!/usr/bin/python

import argparse, sys, os, json, subprocess, re
import datetime
from multiprocessing import cpu_count, Pool
from random import randint
from shutil import copyfile, move, rmtree
from SequenceBasics import GenericFastaFileReader, GenericFastqFileReader
from SequenceCompressionBasics import HomopolymerCompressionFactory
from SamToNavBasics import SamToNavFactory
from NavToMapBasics import NavToMapFactory
from CorrectFromMapBasics import CorrectFromMapFactory
from AlignmentBasics import GenericAlignerCaller, GenericAlignerIndexBuilder

def main():
  # Run modes for runLSC.py are
  # a. Mode 0: Run all the way through (Modes 1-3)
  # b. Mode 1: Run the homopolymer compression to prepare short reads data.
  # c. Mode 2: Run the alignment all the way through correction.
  # c(alternative). Parallelized Mode 2: Run the alignment all the way through correction on a subset of long reads.
  # d. Mode 3: Compile all the individual corrected reads from mode 2 into a single output.
  version = "2.alpha"
  parser = argparse.ArgumentParser(description="LSC "+version+": Correct errors (e.g. homopolymer errors) in long reads, using short read data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--long_reads',required=True,help="FASTAFILE Long reads to correct")
  parser.add_argument('--short_reads',nargs='+',help="FASTA/FASTQ FILE Short reads used to correct the long reads. Can be multiple files.  If choice is cps reads, then there must be 2 files, the cps and the idx file following --short reads")
  parser.add_argument('--short_read_file_type',default='fa',choices=['fa','fq','cps'],help="Short read file type")
  parser.add_argument('--threads',type=int,default=0,help="Number of threads (Default = cpu_count)")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default='/tmp',help="FOLDERNAME where temporary files can be placed")
  group.add_argument('--specific_tempdir',help="FOLDERNAME of exactly where to place temproary folders")
  parser.add_argument('-o','--output',required=True,help="FOLDERNAME where output is to be written")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--mode',default=0,choices=[0,1,2,3],type=int,help="0: run through")
  group.add_argument('--parallelized_mode_2',type=int,help="Mode 2, but you specify a sigle batch to execute.")
  parser.add_argument('--aligner',default='hisat',choices=['hisat','bowtie2'],help="Aligner choice")
  parser.add_argument('--sort_mem_max',type=int,help="-S option for memory in unix sort")
  parser.add_argument('--minNumberofNonN',type=int,default=40,help="Minimum number of non-N characters in the compressed read")
  parser.add_argument('--maxN',type=int,help="Maximum number of Ns in the compressed read")
  parser.add_argument('--error_rate_threshold',type=int,default=12,help="Maximum percent of errors in a read to use the alignment")
  parser.add_argument('--short_read_coverage_threshold',type=int,default=20,help="Minimum short read coverage to do correction")
  parser.add_argument('--long_read_batch_size',type=int,default=5000,help="INT number of long reads to work on at a time")
  parser.add_argument('--samtools_path',default='samtools',help="Path to samtools by default assumes its installed.")
  args = parser.parse_args()
  if args.threads == 0:
    args.threads = cpu_count()

  run_pathfilename = os.path.realpath(__file__)
  bin_path, run_filename = GetPathAndName(run_pathfilename)

  # Establish running parameters here.  Some of these may need to become command line options.
  mode = 0
  if args.parallelized_mode_2:
    mode = 2
  elif args.mode:
    mode = args.mode
  if (mode == 1 or mode == 2) and not args.specific_tempdir:
    sys.stderr.write("ERROR: if you want to run mode 1 or 2, you will need to define a specific temporary directory to store intermediate files with --specific_tempdir FOLDERNAME\n")
    sys.exit()
  LR_filetype = 'fa'
  I_nonredundant = "N"
  sort_max_mem = -1
  clean_up = 0
  bowtie2_options = "--end-to-end -a -f -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.08 --no-unal --omit-sec-seq"

  #Figure out the path for samtools.  If one is not specified, just try the local.
  ofnull = open('/dev/null','w')
  if args.samtools_path == "samtools": #user is not defining their own path
    cmd_sam = 'which samtools'
    p = subprocess.Popen(cmd_sam.split(),stderr=ofnull,stdout=subprocess.PIPE)
    output, err = p.communicate()
    if len(output.rstrip()) > 0 and re.search('samtools$',output.rstrip()):
      #samtools is already installed
      args.samtools_path = output.rstrip()
    else: # try to use a local samtools
      args.samtools_path = bin_path.rstrip('/')+'/../samtools-1.2/samtools'
  ofnull.close()
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

  if mode > 1 and not args.specific_tempdir:
    sys.stderr.write("Error: to run mode 2 or 3 you need to specify a temporary directory with --specific_tempdir.  An easy way to create one is to run mode 1 with the --specific_tempdir option.\n")

  if not os.path.isdir(temp_foldername):
    if mode == 2:
      sys.stderr.write("Error: temp folder does not exist.\n")
      sys.stderr.write("Note: You need to run LSC in mode 1 or 0 before running in mode 2.\n")
      sys.exit()
    else: 
      os.makedirs(temp_foldername)
  t0 = datetime.datetime.now()

  # LSC does have an option for a 'clean_up' step with clean_up.py but I don't recall it being documented 
  # note that clean_up.py took as inputs the temp folder, and the two cpu_counts
  
  # Make sure some folders that will eventually need to be there get there ASAP   
  # These will also get made on the fly if they happen to not be there but we might as well try to make them prior to parallelization
  if not os.path.exists(temp_foldername+'LR_Compressed'):
    os.makedirs(temp_foldername+'LR_Compressed')
  if not os.path.exists(temp_foldername+'Aligner_Indices'):
    os.makedirs(temp_foldername+'Aligner_Indices')
  if not os.path.exists(temp_foldername+'Alignments'):
    os.makedirs(temp_foldername+'Alignments')
  if not os.path.exists(temp_foldername+'Nav_Files'):
    os.makedirs(temp_foldername+'Nav_Files')
  if not os.path.exists(temp_foldername+'Map_Files'):
    os.makedirs(temp_foldername+'Map_Files')
  if not os.path.exists(temp_foldername+'Output_Files'):
    os.makedirs(temp_foldername+'Output_Files')
  if not os.path.exists(temp_foldername+'Log_Files'):
    os.makedirs(temp_foldername+'Log_Files')

  ################################################################################
  # Remove duplicate short reads first  
  if mode == 0 or mode == 1:    
    if I_nonredundant == "N" and args.short_read_file_type != "cps":  # If we go in, we want to get a unique set
        remove_duplicate_short_reads(temp_foldername,args)
        sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")
        
  ############################################################################
  # Run compression over short reads
  if (mode == 0 or mode == 1) and args.short_read_file_type != "cps":
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
        
  #Handle the case where the input is a set of cps/idx files
  #In this case we will have alread skipped compression and making things unique
  if (mode == 0 or mode == 1) and args.short_read_file_type == 'cps':
    if len(args.short_reads) != 2:
      sys.stderr.write("ERROR: Short ready type is cps, but the two inputs, .cps and .idx were not given.\n")
      sys.exit()
    copyfile(args.short_reads[0],temp_foldername+'SR.fa.cps')
    copyfile(args.short_reads[1],temp_foldername+'SR.fa.idx')

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
      of_lr.write(">"+str(z)+"\n"+entry['seq'].upper()+"\n")
      of_lr_readnames.write(str(z)+"\t"+entry['name']+"\n")
    of_lr.close()
    of_lr_readnames.close()
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")

  ######################################################################
  # Find the total number of batches we will run (regardless of run type)
  if (mode == 0 or mode == 1 or mode == 2) and not args.parallelized_mode_2:
    sys.stderr.write("===batch count:===\n")    
    total_batches = batch_count_LR(temp_foldername+'LR.fa.readnames',args.long_read_batch_size)
    of_batch = open(temp_foldername+'batch_count','w')
    of_batch.write(str(total_batches)+"\n")
    of_batch.close()
    sys.stderr.write("Work will begin on "+str(total_batches)+" batches.\n")
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")
  elif mode == 3:
    total_batches = 0
    with open(temp_foldername+"batch_count") as inf:
      read_batches = int(inf.readline().rstrip())
      if not read_batches:
        sys.stderr.write("ERROR. problem reading batch count\n")
        sys.exit()
      total_batches = read_batches
  elif args.parallelized_mode_2:
    total_batches = 0
    with open(temp_foldername+"batch_count") as inf:
      read_batches = int(inf.readline().rstrip())
      if not read_batches:
        sys.stderr.write("ERROR. problem reading batch count\n")
        sys.exit()
      total_batches = read_batches
    if args.parallelized_mode_2 > 0 and args.parallelized_mode_2 <= total_batches:
      sys.stderr.write("===Parallelized mode:===\nWill begin work on batch #"+str(args.parallelized_mode_2)+"\n")
    else:
      sys.stderr.write("ERROR invalid batch number for parallelized_mode_2\n")
      sys.exit()
  sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")

  ##########################################
  # Compress the long reads
  # Build the aligner index 
  # This is the step where parallelization aught to come in to play
  # so I am ending mode 1 at this step.
  if mode == 0 or mode == 2:
    sys.stderr.write("===compress LR, samParser LR.??.cps.nav:===\n")
    if not os.path.exists(temp_foldername+'LR_Compressed'):
      os.makedirs(temp_foldername+'LR_Compressed')
    if not os.path.exists(temp_foldername+'Aligner_Indices'):
      os.makedirs(temp_foldername+'Aligner_Indices')
    execute_pre_alignment(temp_foldername+'LR.fa',args,temp_foldername,args.long_read_batch_size)
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")

  ###############################################
  # Add in a step here for aligning the reads batch by batch
  # one at a time since disk reads on the files seem to come at a premium
  #step 3 map short reads
  if mode == 0 or mode == 2:
    if not os.path.exists(temp_foldername+'/Alignments'):
      os.makedirs(temp_foldername+'/Alignments')
    for batch_number in range(1,total_batches+1):
      if not args.parallelized_mode_2 or args.parallelized_mode_2 == batch_number:
        sys.stderr.write("... step 3 aligning "+str(batch_number)+"/"+str(total_batches)+"   \r")
        of = open(temp_foldername+'/Alignments/LR.fa.'+str(batch_number)+'.cps.bam','w')
        gac = GenericAlignerCaller(args.aligner,temp_foldername+'SR.fa.cps')
        gac.set_samtools_path(args.samtools_path)
        #sys.stderr.write(args.samtools_path+"\n")
        gac.set_threads(args.threads)
        gac.execute(temp_foldername+'Aligner_Indices/'+str(batch_number)+'/LR.fa.'+str(batch_number)+'.cps.aligner-index',of)
        of.close()
    sys.stderr.write("\n")

  ###############################################
  # Add in a step here for returning to the per batch run
  # Make the maps
  # Make the corrections
  if mode == 0 or mode == 2:
    sys.stderr.write("===compress LR, samParser LR.??.cps.nav:===\n")    
    execute_LR(args,temp_foldername,total_batches)
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")
  
  ##########################################
  # Combine the batch outputs into final outputs
  if mode == 0 or mode == 3:
    sys.stderr.write("===produce outputs:===\n")    
    of_uncorrected = open(temp_foldername+"uncorrected_LR.fa",'w')
    of_corrected = open(temp_foldername+"corrected_LR.fa",'w')
    of_corrected_fq = open(temp_foldername+"corrected_LR.fq",'w')
    of_full = open(temp_foldername+"full_LR.fa",'w')
    # read total batches from the file
    if not os.path.isfile(temp_foldername+"batch_count"):
      sys.stderr.write("ERROR you must run in modes 1 and 2 (or 0) before running mode 3\n")
      sys.exit()
    with open(temp_foldername+"batch_count") as inf:
      read_batches = int(inf.readline().rstrip())
      if not read_batches:
        sys.stderr.write("ERROR. problem reading batch count\n")
        sys.exit()
      total_batches = read_batches
    sys.stderr.write("Producing outputs from "+str(total_batches)+" batches\n")
    for i in range(1,total_batches+1):
      if not os.path.isfile(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_uncorrected"):
        sys.stderr.write("Warning missing file "+temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_uncorrected"+"\n")
      else:
        with open(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_uncorrected") as inf:
          for line in inf:
            of_uncorrected.write(line)
      if not os.path.isfile(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_corrected"):
        sys.stderr.write("Warning missing file "+temp_foldername+"output."+str(i)+".file_uncorrected"+"\n")
      else:
        with open(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_corrected") as inf:
          for line in inf:
            of_corrected.write(line)
      if not os.path.isfile(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_corrected_fq"):
        sys.stderr.write("Warning missing file "+temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_corrected_fq"+"\n")
      else:
        with open(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_corrected_fq") as inf:
          for line in inf:
            of_corrected_fq.write(line)
      if not os.path.isfile(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_full"):
        sys.stderr.write("Warning missing file "+temp_foldername+"output."+str(i)+".file_full"+"\n")
      else:
        with open(temp_foldername+'Output_Files/'+str(i)+"/output."+str(i)+".file_full") as inf:
          for line in inf:
            of_full.write(line)
    of_uncorrected.close()
    of_corrected.close()
    of_corrected_fq.close()
    of_full.close()
    copyfile(temp_foldername+"uncorrected_LR.fa",output_foldername+"uncorrected_LR.fa")
    copyfile(temp_foldername+"corrected_LR.fa",output_foldername+"corrected_LR.fa")
    copyfile(temp_foldername+"corrected_LR.fq",output_foldername+"corrected_LR.fq")
    copyfile(temp_foldername+"full_LR.fa",output_foldername+"full_LR.fa")
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")

def GetPathAndName(pathfilename):
  full_path = os.path.abspath(pathfilename)
  [path, filename] = os.path.split(full_path)
  path += '/'
  return path, filename

# Pre: temp_foldername and the input arguments
# Post: Writes SR_uniq.fa into tempfolder as a 

def remove_duplicate_short_reads(temp_foldername,args):
  SR_filetype = args.short_read_file_type
  sys.stderr.write("=== sort and uniq SR data ===\n")
  of = open(temp_foldername+'SR_uniq.seq','w')
  for SR_pathfilename in args.short_reads:
    if SR_filetype == 'fa':
      gfr = GenericFastaFileReader(SR_pathfilename)
    elif SR_filetype == 'fq':
      gfr = GenericFastqFileReader(SR_pathfilename)
    else:
      sys.stderr.write("ERROR: invalid SR_filetype "+SR_filetype+"\n")
      sys.exit()
    #Launch a pipe to store the unique fasta.  Its a little cumbersome but should
    #minimize memory usage.  Could go back and add the -S to sort if memory is a problem
    cmd = 'sort -T '+temp_foldername
    if args.sort_mem_max:
      cmd += " -S "+str(args.sort_mem_max)
    cmd += " | uniq -c "
    p = subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=of,shell=True)
    while True:
      entry = gfr.read_entry()
      if not entry: break
      seq = entry['seq'].upper()
      p.stdin.write(seq+"\n")
    p.communicate()
  of.close()
  z = 0
  output_SR_pathfilename = temp_foldername + "SR_uniq.fa"
  of = open(output_SR_pathfilename,'w')
  with open(temp_foldername+'SR_uniq.seq') as inf:
    for line in inf:
      z += 1
      line = line.rstrip()
      m = re.match('^\s*(\d+)\s+(\S+)$',line)
      if not m:
        sys.stderr.write("ERROR unable to process uniq -c line "+line+"\n")
      of.write(">"+str(z)+"_"+m.group(1)+"\n"+m.group(2)+"\n")
  of.close
  return output_SR_pathfilename

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

#########################################
# see based on our 
def batch_count_LR(lr_filename,batchsize):
  gfr = GenericFastaFileReader(lr_filename)
  batch_current = 0
  batch_number = 0
  with open(lr_filename) as inf:
    for line in inf:
      batch_current += 1
      if batch_current >= batchsize:
        batch_number += 1
        batch_current = 0
    if batch_current > 0:
      batch_number += 1
  return batch_number

# This function launches our jobs (mode 2), and decides at which level (if any) to parallelize the job based on mode and thread counts.
# Pre: long read files
#      input arguments
#      temporary folder name
#      batch size
# Post: For each batch of batch size that runs, the temporarly folder
#       will get the results of LSC for that batch.  
#       Files output are 
#     output.#.file_uncorrected
#     output.#.file_uncorrected_fq
#     output.#.file_corrected
#     output.#.file_full
#       where # is the batch number
# Modifies: Writes temporary files. Its tricky in how it runs.
#   If you're using a single computer to run through and its not in the --parallelized_mode_2,
#     it will run a Pool of jobs and give a single core to each job to run over the batch_size of long reads
#   If you're using a single computer to run a parallelized_mode_2 job where you going to execute one of 
#     the sets from the batches, it will use the number of threads you specify within that batch to speed processing
#   The steps it runs through are better defined in execute_batch 
def execute_LR(args,temp_foldername,total_batches):
  if args.threads > 1 and not args.parallelized_mode_2: 
    p = Pool(processes=args.threads)
  for batch_number in range(1,total_batches+1):
    # We will not use Pool if we are only using one thread or we are doing parallelized mode
    if (not args.parallelized_mode_2 and args.threads <= 1) or args.parallelized_mode_2 == batch_number:
      execute_batch(batch_number,args.minNumberofNonN,args.maxN,temp_foldername,args.threads,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max,args.samtools_path)
    # We will launch by Pool if we are only running on one system.
    # In this case we only give each pool job one processor
    elif not args.parallelized_mode_2:
      p.apply_async(execute_batch,args=(batch_number,args.minNumberofNonN,args.maxN,temp_foldername,1,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max,args.samtools_path))
  #If we are running jobs as pools we clean up those threads here.
  if args.threads > 1 and not args.parallelized_mode_2:
    p.close()
    p.join()
  sys.stderr.write("\n")
  return

def execute_pre_alignment(lr_filename,args,temp_foldername,batchsize):
  gfr = GenericFastaFileReader(lr_filename)
  batch = []
  batch_number = 0
  if args.threads > 1 and not args.parallelized_mode_2: 
    p = Pool(processes=args.threads)
  while True:
    entry = gfr.read_entry()
    if not entry: break
    batch.append([entry['name'],entry['seq']])
    if len(batch) >= batchsize:
      batch_number += 1
      # We will not use Pool if we are only using one thread or we are doing parallelized mode
      if (not args.parallelized_mode_2 and args.threads <= 1) or args.parallelized_mode_2 == batch_number:
        execute_batch_pre_alignment(json.dumps(batch),batch_number,args.minNumberofNonN,args.maxN,temp_foldername,args.threads,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max,args.aligner)
      # We will launch by Pool if we are only running on one system.
      # In this case we only give each pool job one processor
      elif not args.parallelized_mode_2:
        p.apply_async(execute_batch_pre_alignment,args=(json.dumps(batch),batch_number,args.minNumberofNonN,args.maxN,temp_foldername,1,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max,args.aligner))
      batch = []
  #Handle the case were a buffered sequence(s) still needs run.
  if len(batch) > 0:
    batch_number += 1
    # We will not use Pool if we are only using one thread or we are doing parallelized mode
    if (not args.parallelized_mode_2 and args.threads <= 1) or args.parallelized_mode_2 == batch_number:
      execute_batch_pre_alignment(json.dumps(batch),batch_number,args.minNumberofNonN,args.maxN,temp_foldername,args.threads,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max,args.aligner)
    # We will launch by Pool if we are only running on one system.
    # In this case we will only give each pool job one processor
    elif not args.parallelized_mode_2:
      p.apply_async(execute_batch_pre_alignment,args=(json.dumps(batch),batch_number,args.minNumberofNonN,args.maxN,temp_foldername,1,args.error_rate_threshold,args.short_read_coverage_threshold,args.sort_mem_max,args.aligner))
  #If we are running jobs as pools we clean up those threads here.
  if args.threads > 1 and not args.parallelized_mode_2:
    p.close()
    p.join()
  sys.stderr.write("\n")
  return batch_number

def execute_batch(batch_number,minNumberofNonN,maxN,temp_foldername,threads,error_rate_threshold,short_read_coverage_threshold,sort_max_mem,samtools_path):
  if not os.path.exists(temp_foldername+'Log_Files'):
    os.makedirs(temp_foldername+'Log_Files')
  of_log = open(temp_foldername+"Log_Files/LR.fa."+str(batch_number)+'.cps'+'.log3','w')
  # At this step we can remove the index
  rmtree(temp_foldername+'Aligner_Indices/'+str(batch_number)+'/')

  #step 4 convert sam to nav
  if not os.path.exists(temp_foldername+'Nav_Files'):
    os.makedirs(temp_foldername+'Nav_Files')
  sys.stderr.write("... executing batch "+str(batch_number)+".4 (bam to nav)   \r")
  of = open(temp_foldername+'Nav_Files/LR.fa.'+str(batch_number)+'.cps.nav','w')
  stnf = SamToNavFactory()
  stnf.set_error_rate_threshold(error_rate_threshold)
  stnf.initialize_compressed_query(temp_foldername+'SR.fa.cps',temp_foldername+'SR.fa.idx')
  stnf.initialize_target(temp_foldername+'LR_Compressed/'+str(batch_number)+'/LR.fa.'+str(batch_number)+'.cps')
  cmd4 = samtools_path+' view '+temp_foldername+'/Alignments/LR.fa.'+str(batch_number)+'.cps.bam'
  p = subprocess.Popen(cmd4,stdout=subprocess.PIPE,shell=True)
  while True:
    line = p.stdout.readline()
    if not line: break
    navs = stnf.sam_to_nav(line)
    if not navs: continue
    for oline in navs:
      of.write(oline+"\n")
  p.communicate()
  stnf.close()
  of.close()

  # At this step we can remove the bam
  os.remove(temp_foldername+"Alignments/LR.fa."+str(batch_number)+".cps.bam")
  
  #step 5 sort nav by long read name
  sys.stderr.write("...executing batch "+str(batch_number)+".5 (sort nav)    \r")
  sort_cmd =  "sort -T "+temp_foldername+" "
  if sort_max_mem:
    sort_cmd += "-S "+str(sort_max_mem)+" "
  sort_cmd += "-nk 2 " + temp_foldername+ "Nav_Files/LR.fa."+str(batch_number)+".cps.nav "
  sort_cmd += "-o "+temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav.sort"
  subprocess.call(sort_cmd.split(),stderr=of_log,stdout=of_log)  

  #step 6 nav to mapping
  if not os.path.exists(temp_foldername+'Map_Files'):
    os.makedirs(temp_foldername+'Map_Files')
  sys.stderr.write("... executing batch "+str(batch_number)+".6 (nav to map)   \r")
  ntmf = NavToMapFactory(temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav.sort", \
         temp_foldername+"LR_Compressed/"+str(batch_number)+"/LR.fa."+str(batch_number)+".cps", \
         temp_foldername+"LR_Compressed/"+str(batch_number)+"/LR.fa."+str(batch_number)+".idx", \
         temp_foldername+"LR.fa.readnames", \
         short_read_coverage_threshold)
  of_map = open(temp_foldername+"Map_Files/LR_SR.map."+str(batch_number),'w')
  while True:
    entry = ntmf.read_entry()
    if not entry: break
    of_map.write(entry+"\n")
  ntmf.close()
  of_map.close()

  # At this step we can remove the nav
  os.remove(temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav")
  os.remove(temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav.sort")

  #step 7 correct
  sys.stderr.write("... executing batch "+str(batch_number)+".7 (correcting)   \r")
  if not os.path.exists(temp_foldername+'Output_Files'):
    os.makedirs(temp_foldername+'Output_Files')
  if not os.path.exists(temp_foldername+'Output_Files/'+str(batch_number)):
    os.makedirs(temp_foldername+'Output_Files/'+str(batch_number))
  output_prefix = temp_foldername+'Output_Files/'+str(batch_number)+"/output."+str(batch_number)+".file"
  temp_prefix = temp_foldername+'Output_Files/'+str(batch_number)+"/tempoutput."+str(batch_number)+".file"
  full_read_file=open(temp_prefix+'_full','w')
  corrected_read_file=open(temp_prefix+ '_corrected','w')
  corrected_read_fq_file=open(temp_prefix+'_corrected_fq','w')
  uncorrected_read_file = open(temp_prefix+'_uncorrected','w')
  cfmf = CorrectFromMapFactory(temp_foldername+"LR.fa.readnames")
  with open(temp_foldername+"Map_Files/LR_SR.map."+str(batch_number)) as inf:
    for line in inf:
      output = cfmf.correct_map_line(line)
      if not output: continue
      uncorrected_read_file.write(output['uncorrected_fasta'])
      corrected_read_file.write(output['corrected_fasta'])
      corrected_read_fq_file.write(output['corrected_fastq'])
      full_read_file.write(output['full_fasta'])
  corrected_read_file.close()
  corrected_read_fq_file.close()
  uncorrected_read_file.close()
  full_read_file.close()
  move(temp_prefix+'_full',output_prefix+'_full')
  move(temp_prefix+'_uncorrected',output_prefix+'_uncorrected')
  move(temp_prefix+'_corrected',output_prefix+'_corrected')
  move(temp_prefix+'_corrected_fq',output_prefix+'_corrected_fq')

  of_log.close()
  return

def execute_batch_pre_alignment(batch_json,batch_number,minNumberofNonN,maxN,temp_foldername,threads,error_rate_threshold,short_read_coverage_threshold,sort_max_mem,aligner):
  batch = json.loads(batch_json)
  #step 1 compression
  sys.stderr.write("... executing batch "+str(batch_number)+".1 (compression)   \r")
  if not os.path.exists(temp_foldername+'LR_Compressed/'+str(batch_number)):
    os.makedirs(temp_foldername+'LR_Compressed/'+str(batch_number))
  hpcf = HomopolymerCompressionFactory()
  if minNumberofNonN: hpcf.set_MinNonN(minNumberofNonN)
  if maxN: hpcf.set_MaxN(maxN)
  of_cps = open(temp_foldername+"LR_Compressed/"+str(batch_number)+"/LR.fa."+str(batch_number)+'.cps','w')
  of_idx = open(temp_foldername+"LR_Compressed/"+str(batch_number)+"/LR.fa."+str(batch_number)+'.idx','w')
  for entry in batch:
    res = hpcf.compress(entry[1])
    if not res: continue
    (cps, pos_ls, len_ls) = res
    of_cps.write(">"+entry[0]+"\n"+cps+"\n")
    of_idx.write(entry[0]+"\t"+pos_ls+"\t"+len_ls+"\n")
  of_cps.close()
  of_idx.close()

  #step 2 build an index
  sys.stderr.write("... executing batch "+str(batch_number)+".2 (index)          \r")
  if not os.path.exists(temp_foldername+'Aligner_Indices/'+str(batch_number)):
    os.makedirs(temp_foldername+'Aligner_Indices/'+str(batch_number))
  infile = temp_foldername+'LR_Compressed/'+str(batch_number)+'/LR.fa.'+str(batch_number)+'.cps'
  gaib = GenericAlignerIndexBuilder(aligner,infile)
  gaib.execute(temp_foldername+'Aligner_Indices/'+str(batch_number)+'/LR.fa.'+str(batch_number)+'.cps.aligner-index')
  return

if __name__=="__main__":
  main()
