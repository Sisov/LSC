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
  version = "2.0"
  parser = argparse.ArgumentParser(description="LSC "+version+": Correct errors (e.g. homopolymer errors) in long reads, using short read data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--long_reads',help="FASTAFILE Long reads to correct. Required in mode 0 or 1.")
  parser.add_argument('--short_reads',nargs='*',help="FASTA/FASTQ FILE Short reads used to correct the long reads. Can be multiple files.  If choice is cps reads, then there must be 2 files, the cps and the idx file following --short reads. Required in mode 0 or 1.")
  parser.add_argument('--short_read_file_type',default='fa',choices=['fa','fq','cps'],help="Short read file type")
  parser.add_argument('--threads',type=int,default=0,help="Number of threads (Default = cpu_count)")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default='/tmp',help="FOLDERNAME where temporary files can be placed")
  group.add_argument('--specific_tempdir',help="FOLDERNAME of exactly where to place temproary folders. Required in mode 1, 2 or 3. Recommended for any run where you may want to look back at intermediate files.")
  parser.add_argument('-o','--output',help="FOLDERNAME where output is to be written. Required in mode 0 or 3.")
  group1 = parser.add_mutually_exclusive_group()
  group1.add_argument('--mode',default=0,choices=[0,1,2,3],type=int,help="0: run through, 1: Prepare homopolymer compressed long and short reads.  2: Execute correction on batches of long reads.  Can be superseded by --parallelized_mode_2 where you will only execute a single batch.  3: Combine corrected batches into a final output folder.")
  group1.add_argument('--parallelized_mode_2',type=int,help="Mode 2, but you specify a sigle batch to execute.")
  parser.add_argument('--aligner',default='bowtie2',choices=['hisat','bowtie2'],help="Aligner choice. hisat parameters have not been optimized, so we recommend bowtie2.")
  parser.add_argument('--sort_mem_max',type=int,help="-S option for memory in unix sort")
  parser.add_argument('--minNumberofNonN',type=int,default=40,help="Minimum number of non-N characters in the compressed read")
  parser.add_argument('--maxN',type=int,help="Maximum number of Ns in the compressed read")
  parser.add_argument('--error_rate_threshold',type=int,default=12,help="Maximum percent of errors in a read to use the alignment")
  parser.add_argument('--short_read_coverage_threshold',type=int,default=20,help="Minimum short read coverage to do correction")
  parser.add_argument('--long_read_batch_size',default=500,type=int,help="INT number of long reads to work on at a time.  This is a key parameter to adjusting performance.  A smaller batch size keeps the sizes and runtimes of intermediate steps tractable on large datasets, but can slow down execution on small datasets.  The default value should be suitable for large datasets.")
  parser.add_argument('--samtools_path',default='samtools',help="Path to samtools by default assumes its installed.  If not specified, the included version will be used.")
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

  check_argument_mode_compatibility(mode,args)

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
  if not os.path.exists(temp_foldername+'Split_Short'):
    os.makedirs(temp_foldername+'Split_Short')
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

  ################################################################################
  # Remove duplicate short reads first  
  if mode == 0 or mode == 1:    
    if args.short_read_file_type != "cps":  # If we go in, we want to get a unique set
      remove_duplicate_short_reads(temp_foldername,args)
      sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")

  ############################################################################
  # Run compression over short reads
  if (mode == 0 or mode == 1) and args.short_read_file_type != "cps":
    sys.stderr.write("===compress SR:===\n")    
    # At this point we should have uniq fasta formatted reads
    # We may want to consider supportin a unique set of fastq reads also
    compress_SR(temp_foldername,args)
    sys.stderr.write(str(datetime.datetime.now()-t0)+"\n")
        
  #Handle the case where the input is a set of cps/idx files
  #In this case we will have alread skipped compression and making things unique
  if (mode == 0 or mode == 1) and args.short_read_file_type == 'cps':
    if len(args.short_reads) != 2:
      sys.stderr.write("ERROR: Short ready type is cps, but the two inputs, .cps and .idx were not given.\n")
      sys.exit()
    copyfile(args.short_reads[0],temp_foldername+'SR.fa.cps')
    copyfile(args.short_reads[1],temp_foldername+'SR.fa.idx')

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
    output_foldername = args.output.rstrip('/')+'/'
    if not os.path.isdir(output_foldername):
      os.makedirs(output_foldername)
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

#Pre:  A path to a file
#Post: Return the directory path (with a trailing forward slash), and filename
def GetPathAndName(pathfilename):
  full_path = os.path.abspath(pathfilename)
  [path, filename] = os.path.split(full_path)
  path += '/'
  return path, filename

# Pre: temp_foldername and the input arguments
# Post: Writes SR_uniq.fa into tempfolder as fasta with headers in the format to SR_uniq.fa in the temp directory
#       >X_Y  Where X is a unique integer and Y is the number of times the sequence was observed in the original data
# Modifies: Makes a temporary SR_uniq.seq that has only the unique sequences (sorted) along the way.
def remove_duplicate_short_reads(temp_foldername,args):
  SR_filetype = args.short_read_file_type
  sys.stderr.write("=== sort and uniq SR data ===\n")
  cwd = os.path.dirname(os.path.realpath(__file__))
  udir = cwd+'/../utilities'

  cmd4 = 'split -l 5000000 - '+temp_foldername+'Split_Short/unsort.'
  p4 = subprocess.Popen(cmd4.split(),stdin=subprocess.PIPE,bufsize=1)
  for SR_pathfilename in args.short_reads:
    #Running these through command line scripts because its faster than python
    cmd1 = 'cat '+SR_pathfilename
    if re.search('\.gz$',SR_pathfilename):
      cmd1 = 'zcat '+SR_pathfilename
    if SR_filetype == 'fa':
      cmd2 = udir+'/fasta_to_tsv.pl'
    elif SR_filetype == 'fq':
      cmd2 = udir+'/fastq_to_tsv.pl'
    else:
      sys.stderr.write("ERROR: invalid SR_filetype "+SR_filetype+"\n")
      sys.exit()
    cmd3 = 'cut -f 2'
    p3 = subprocess.Popen(cmd3.split(),stdin=subprocess.PIPE,stdout=p4.stdin,bufsize=1)
    p2 = subprocess.Popen(cmd2.split(),stdin=subprocess.PIPE,stdout=p3.stdin,bufsize=1)
    p1 = subprocess.Popen(cmd1.split(),stdout=p2.stdin,bufsize=1)
    p1.communicate()
    p2.communicate()
    p3.communicate()
    # reads through all SR files this way
  p4.communicate()

  # Now we do the sort operation
  of = open(temp_foldername+'SR_uniq.seq','w')
  fnames = os.listdir(temp_foldername+'Split_Short')
  if args.threads > 1:
    p = Pool(processes=args.threads)
  snames = []
  for n in fnames:
    path = temp_foldername+'Split_Short/'+n
    m = re.search('^(.*)unsort\.([^\.]+)$',path)
    if not m: continue # must be a left over sorted output from another run
    spath = m.group(1)+'sort.'+m.group(2)
    snames.append(spath)
    if args.threads > 1:
      p.apply_async(do_sort,args=(path,spath,))
    else:
      do_sort(path,spath)
    #print path
    #print spath
  if args.threads > 1:
    p.close()
    p.join()
  of.close()

  #### Temporary can delete when not debuging (should be fine even if you dont)
  #snames = [temp_foldername+'Split_Short/'+x for x in os.listdir(temp_foldername+'Split_Short') if re.match('^sort',x)]

  # now we can put them back together
  if len(snames)==0:
    sys.stderr.write("ERROR no sorted short reads\n")
    sys.exit()
  # now final merge
  cmd1 = 'sort --batch-size='+str(max(2,len(snames)))+' -m'
  for n in snames:
    cmd1 += ' '+n
  of = open(temp_foldername+'SR_uniq.seq','w')
  p1 = subprocess.Popen(cmd1.split(),stdout=subprocess.PIPE)
  cmd2 = 'uniq -c'
  p2 = subprocess.Popen(cmd2.split(),stdin=p1.stdout,stdout=of)
  p2.communicate()
  p1.communicate()
  of.close()
  for n in snames: os.remove(n)
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


#Simple call to sort
def do_sort(path,spath):
  dn = open(os.devnull,'w')
  of = open(spath,'w')
  cmd1 = 'sort '+path
  p = subprocess.Popen(cmd1.split(),stdout=of,bufsize=1,stderr=dn)
  p.communicate()
  os.remove(path)
  return

# Handle homopolymer compressing all the short read sequences
# will work on compressing sequences in batches of the batchsize defined in the function
# Pre: The short read file (SR_uniq.fa), and the input arguments
#       globals SR_cps_fh and SR_idx_fh need to point to writable files
# Post: Writes out the individusl results
# Modifies: Launch jobs in parallel through Pool
def compress_SR(temp_filename,args):
  #Lets split the reads to compress
  short_read_file = temp_filename+'SR_uniq.fa'
  cwd = os.path.dirname(os.path.realpath(__file__))
  udir = cwd+'/../utilities'
  fcount = None

  cmd = udir+'/explode_fasta.pl 1000000 '+temp_filename+'/Split_Short'
  inf = open(short_read_file)
  p = subprocess.Popen(cmd.split(),stdin=inf,stdout=subprocess.PIPE)
  fcount = int(p.communicate()[0].rstrip())

  if not fcount:
    fcount = len([x for x in os.listdir(temp_filename+'/Split_Short') if re.match('^\d+\.fa$',x)])

  #sys.stderr.write(str(fcount)+" blocks\n")
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for i in range(1,fcount+1):
    if args.threads > 1:
      p.apply_async(compress_batch,args=(temp_filename,i,args.minNumberofNonN,args.maxN,))
    else:
      compress_batch(temp_filename,i,args.minNumberofNonN,args.maxN)

  if args.threads > 1:
    p.close()
    p.join()

  #rejoin files
  SR_cps_fh = open(temp_filename+'SR.fa.cps','w')
  SR_idx_fh = open(temp_filename+'SR.fa.idx','w')
  for i in range(1,fcount+1):
    fname1 = temp_filename+'Split_Short/'+str(i)+'.fa.cps'
    with open(fname1) as inf:
      for line in inf:  SR_cps_fh.write(line)
    os.remove(fname1)
    fname2 = temp_filename+'Split_Short/'+str(i)+'.fa.idx'
    with open(fname2) as inf:
      for line in inf:  SR_idx_fh.write(line)
    os.remove(fname2)
  SR_cps_fh.close()
  SR_idx_fh.close()

# Compress batch
# Pre: temp_foldername and index of the split file to process, and parameters minNumofNonN and maxN
# Post: name, compressed_sequence, possition array and length array
def compress_batch(temp_foldername,i,minNumberofNonN,maxN):
  fname = temp_foldername+'Split_Short/'+str(i)+'.fa'
  gfr = GenericFastaFileReader(fname)
  seq_batch = []
  hpcf = HomopolymerCompressionFactory()
  if minNumberofNonN: hpcf.set_MinNonN(minNumberofNonN)
  if maxN: hpcf.set_MaxN(maxN)
  of1 = open(temp_foldername+'Split_Short/'+str(i)+'.fa.cps','w')
  of2 = open(temp_foldername+'Split_Short/'+str(i)+'.fa.idx','w')
  while True:
    entry = gfr.read_entry()
    if not entry: break
    res = hpcf.compress(entry['seq'])
    if not res: continue
    (cps, pos_ls, len_ls) = res
    result = [entry['name'],cps,pos_ls,len_ls]
    of1.write(">"+result[0]+"\n"+result[1]+"\n")
    of2.write(result[0]+"\t"+result[2]+"\t"+result[3]+"\n")
  gfr.close()
  of1.close()
  of2.close()
  os.remove(fname)

#########################################
# Get the total number of Batches
# Pre: lr_filename is long read file
#      batchsize is the number of batches
# Post: The number of batches that we have to run given that batch size and number of long reads
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
      execute_batch(batch_number,temp_foldername,args.threads,args,args.samtools_path,args.long_read_batch_size)
    # We will launch by Pool if we are only running on one system.
    # In this case we only give each pool job one processor
    elif not args.parallelized_mode_2:
      p.apply_async(execute_batch,args=(batch_number,temp_foldername,1,args,args.samtools_path,args.long_read_batch_size))
  #If we are running jobs as pools we clean up those threads here.
  if args.threads > 1 and not args.parallelized_mode_2:
    p.close()
    p.join()
  sys.stderr.write("\n")
  return

# Launch parallel jobs for pre-alignment steps.  These will work on batches of long reads
# Pre: Long read file, arugments, temporary folder, batchsize
# Post: Writes to the temporary file compressed long reads (for each batch), 
#       and an aligner index file based on each batch of long reads
# Modifies: Parallel version of execute_batch_pre_alignment are launched.
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

#For a given batch number write the nav file based on the alignment bam file and alignments
# Pre: temporary folder, batch number, job sub-index, error_rate_threshold, samtools_path
# Post: Nav file for the batch number and sub-index is created.
def write_nav(temp_foldername,batch_number,threads,i,error_rate_threshold,samtools_path):
  of = open(temp_foldername+'Nav_Files/LR.fa.'+str(batch_number)+'.cps.'+str(i)+'.nav','w')
  stnf = SamToNavFactory()
  stnf.set_thread_count(threads)
  stnf.set_thread_index(i)
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

# execute_batch works with a set of long reads to that has already 
# been prepared into alignment indexs and is accessed according to batch_number
# Pre: Requires an alignment index already be present in the temporary folder for each batch from 1 to batch_number
#      batch_number tells it which batch number this function will be working on
#      temp_foldername is the temporary folder we are working in
#      error_rate_threshold is an acceptable error rate
#      short_read_coverage_threshold is an acceptable short read coverage
#      sort_max_mem is the maximum memory to be used during the sort by a thread
#      samtools_path is the location of samtools
#      long_read_batch_size is the number of long reads in each batch
# Post: Creats nav and map files along the way but ultimately stores corrected outputs in the temporary folder 
# Modifies: Temp files, can execute steps in parallel if parallelized_mode_2 is being run and multiple threads are specified.
def execute_batch(batch_number,temp_foldername,threads,args,samtools_path,long_read_batch_size):
  if not os.path.exists(temp_foldername+'Log_Files'):
    os.makedirs(temp_foldername+'Log_Files')
  of_log = open(temp_foldername+"Log_Files/LR.fa."+str(batch_number)+'.cps'+'.log3','w')
  # At this step we can remove the index
  #rmtree(temp_foldername+'Aligner_Indices/'+str(batch_number)+'/')

  #step 4 convert sam to nav
  if not os.path.exists(temp_foldername+'Nav_Files'):
    os.makedirs(temp_foldername+'Nav_Files')
  sys.stderr.write("... executing batch "+str(batch_number)+".4 (bam to nav)   \r")
  if threads > 1:
    p = Pool(processes = threads)
  for i in range(1,threads+1):
    if threads > 1:
       p.apply_async(write_nav,args=(temp_foldername,batch_number,threads,i,args.error_rate_threshold,samtools_path))
    else:
      write_nav(temp_foldername,batch_number,threads,i,args.error_rate_threshold,samtools_path)
  if threads > 1:
    p.close()
    p.join()

  of = open(temp_foldername+'Nav_Files/LR.fa.'+str(batch_number)+'.cps.nav','w')
  for i in range(1,threads+1):
    with open(temp_foldername+'Nav_Files/LR.fa.'+str(batch_number)+'.cps.'+str(i)+'.nav') as inf:
      for line in inf:
        of.write(line)
    os.remove(temp_foldername+'Nav_Files/LR.fa.'+str(batch_number)+'.cps.'+str(i)+'.nav')
  of.close()

  # At this step we can remove the bam
  #os.remove(temp_foldername+"Alignments/LR.fa."+str(batch_number)+".cps.bam")
  
  #step 5 sort nav by long read name
  sys.stderr.write("...executing batch "+str(batch_number)+".5 (sort nav)    \r")
  sort_cmd =  "sort -T "+temp_foldername+" "
  sort_max_mem = args.sort_mem_max
  if sort_max_mem:
    sort_cmd += "-S "+str(sort_max_mem)+" "
  sort_cmd += "-nk 2 " + temp_foldername+ "Nav_Files/LR.fa."+str(batch_number)+".cps.nav "
  sort_cmd += "-o "+temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav.sort"
  subprocess.call(sort_cmd.split(),stderr=of_log,stdout=of_log)  
  # At this step we can remove the nav
  os.remove(temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav")

  #step 6 nav to mapping
  if not os.path.exists(temp_foldername+'Map_Files'):
    os.makedirs(temp_foldername+'Map_Files')
  sys.stderr.write("... executing batch "+str(batch_number)+".6 (nav to map)   \r")
  ntmf = NavToMapFactory(temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav.sort", \
         temp_foldername+"LR_Compressed/"+str(batch_number)+"/LR.fa."+str(batch_number)+".cps", \
         temp_foldername+"LR_Compressed/"+str(batch_number)+"/LR.fa."+str(batch_number)+".idx", \
         temp_foldername+"LR.fa.readnames", \
         args.short_read_coverage_threshold)
  of_map = open(temp_foldername+"Map_Files/LR_SR.map."+str(batch_number),'w')
  while True:
    entry = ntmf.read_entry()
    if not entry: break
    of_map.write(entry+"\n")
  ntmf.close()
  of_map.close()

  # At this step we can remove the nav
  os.remove(temp_foldername+"Nav_Files/LR.fa."+str(batch_number)+".cps.nav.sort")

  #step 7 correct
  sys.stderr.write("... executing batch "+str(batch_number)+".7 (correcting)   \r")
  if not os.path.exists(temp_foldername+'Output_Files'):
    os.makedirs(temp_foldername+'Output_Files')
  if not os.path.exists(temp_foldername+'Output_Files/'+str(batch_number)):
    os.makedirs(temp_foldername+'Output_Files/'+str(batch_number))
  output_prefix = temp_foldername+'Output_Files/'+str(batch_number)+"/output."+str(batch_number)+".file"
  temp_prefix = temp_foldername+'Output_Files/'+str(batch_number)+"/tempoutput."+str(batch_number)+".file"

  buffer = []
  if threads > 1:
    p = Pool(processes=threads)
  correction_batch_number = 0
  long_read_batch_size = int(long_read_batch_size/(threads*100))+1 # work on a few reads at a time
  with open(temp_foldername+"Map_Files/LR_SR.map."+str(batch_number)) as inf:
    for line in inf:
      buffer.append(line)
      if len(buffer) >= long_read_batch_size:
        correction_batch_number += 1
        if threads > 1:
          p.apply_async(do_corrections,args=(buffer,temp_foldername,correction_batch_number,temp_prefix))
        else:
          do_corrections(buffer,temp_foldername,correction_batch_number,temp_prefix)
        buffer = []
  if len(buffer) > 0:
    correction_batch_number += 1
    if threads > 1:
      p.apply_async(do_corrections,args=(buffer,temp_foldername,correction_batch_number,temp_prefix))
    else:
      do_corrections(buffer,temp_foldername,correction_batch_number,temp_prefix)
  if threads > 1:
    p.close()
    p.join()
  of1 = open(temp_prefix+'_full','w')
  of2 = open(temp_prefix+'_uncorrected','w')
  of3 = open(temp_prefix+'_corrected','w')
  of4 = open(temp_prefix+'_corrected_fq','w')
  for num in range(1,correction_batch_number+1):
    with open(temp_prefix+'_full.'+str(num)) as inf:
      for line in inf:
        of1.write(line)
    with open(temp_prefix+'_uncorrected.'+str(num)) as inf:
      for line in inf:
        of2.write(line)
    with open(temp_prefix+'_corrected.'+str(num)) as inf:
      for line in inf:
        of3.write(line)
    with open(temp_prefix+'_corrected_fq.'+str(num)) as inf:
      for line in inf:
        of4.write(line)
  of1.close()
  of2.close()
  of3.close()
  of4.close()
  move(temp_prefix+'_full',output_prefix+'_full')
  move(temp_prefix+'_uncorrected',output_prefix+'_uncorrected')
  move(temp_prefix+'_corrected',output_prefix+'_corrected')
  move(temp_prefix+'_corrected_fq',output_prefix+'_corrected_fq')

  of_log.close()
  return

# Execute corrections on a batch of reads
# Pre: batch - is an array of map file lines ready for correction
#      temp_foldername the temporary foldeer we are working in
#      correction_batch_number the subset of the job we are working on 
#      temp_prefix the file prefix for where we are storing our temporary files
# Post: Written to temp_prefix we have full, corrected, uncorrected fasta/fastq files.
# Modifies: Temporary folder contents
def do_corrections(batch,temp_foldername,correction_batch_number,temp_prefix):
  cfmf = CorrectFromMapFactory(temp_foldername+"LR.fa.readnames")
  temp = {}
  full_read_file=open(temp_prefix+'_full.'+str(correction_batch_number),'w')
  corrected_read_file=open(temp_prefix+ '_corrected.'+str(correction_batch_number),'w')
  corrected_read_fq_file=open(temp_prefix+'_corrected_fq.'+str(correction_batch_number),'w')
  uncorrected_read_file = open(temp_prefix+'_uncorrected.'+str(correction_batch_number),'w')
  outputs = []
  for line in batch:
    output = cfmf.correct_map_line(line)
    uncorrected_read_file.write(output['uncorrected_fasta'])
    corrected_read_file.write(output['corrected_fasta'])
    corrected_read_fq_file.write(output['corrected_fastq'])
    full_read_file.write(output['full_fasta'])
    if not output: continue
    outputs.append(output)
  full_read_file.close()
  corrected_read_file.close()
  corrected_read_fq_file.close()
  uncorrected_read_file.close()

# Execute the prealignment steps of long read compression then index building
# Pre: batch_json - json enocded batch of long reads
#      minNumberofNonN - smallest number of on N values to consider
#      maxN - largest number of N values to consider.
#      temp_foldername the folder we are working in
# 
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

# Pre our run mode 0, 1, 2 or 3, and our input args from argparse
# Post: Returns nothing. Exits hard with error message if incompatabilities in mode and parameter requirements are present
def check_argument_mode_compatibility(mode,args):
  # If we are mode 0 or mode 1 make sure we have required parameters
  if mode == 0 or mode == 1:
    if not args.long_reads:
      sys.stderr.write("ERROR: please specify --long_reads when in mode 0 or mode 1\n")
      sys.exit()
    if not args.short_reads:
      sys.stderr.write("ERROR: please specify --short_reads when in mode 0 or mode 1\n")
      sys.exit()

  if (mode == 1 or mode == 2) and not args.specific_tempdir:
    sys.stderr.write("ERROR: if you want to run mode 1 or 2, you will need to define a specific temporary directory to store intermediate files with --specific_tempdir FOLDERNAME\n")
    sys.exit()

  if mode == 0 or mode == 3:
    if not args.output:
      sys.stderr.write("ERROR please specify output folder when running mode 0 or 3\n")
      sys.exit()

  if mode > 1 and not args.specific_tempdir:
    sys.stderr.write("Error: to run mode 2 or 3 you need to specify a temporary directory with --specific_tempdir.  An easy way to create one is to run mode 1 with the --specific_tempdir option.\n")
    sys.exit()
  return

if __name__=="__main__":
  main()
