import sys, subprocess, os
from multiprocessing import cpu_count

class GenericAlignerCaller:
  def __init__(self,aligner,input_file):
    self.aligner_choice = aligner
    self.output_type = 'bam'
    self.samtools_path = 'samtools' #default is that samtools is installed
    self.input_file = input_file
    self.input_file_type = 'fa'
    self.error_output_handle = open('/dev/null','w')
    self.threads = cpu_count()
    return

  def set_aligner_choice(self,aligner):
    self.aligner_choice = aligner

  def set_threads(self,threads):
    self.threads = threads
    return

  def execute(self,index_file,output_handle):
    self.output_handle = output_handle
    self.index_file = index_file
    if self.aligner_choice == 'hisat':
      self.execute_hisat()
    elif self.aligner_choice == 'bowtie2':
      self.execute_bowtie2()

  def execute_hisat(self):
    cmd3 = ''
    if self.input_file_type == 'fa': cmd3 = ' -f '
    if self.input_file_type == 'fq': cmd3 = ' -q '
    cmd = 'hisat --end-to-end --no-spliced-alignment --no-unal --omit-sec-seq -a --reorder '+cmd3+' -p '+str(self.threads)+' -x '+self.index_file+' -U '+self.input_file
    cmd2 = ' | '+self.samtools_path+' view -Sb -'
    if self.output_type == 'bam':
      cmd += cmd2
    #sys.stderr.write(cmd+"\n")
    subprocess.call(cmd,shell=True,stdout=self.output_handle,stderr=self.error_output_handle)
    return

  def execute_bowtie2(self):
    cmd3 = ''
    if self.input_file_type == 'fa': cmd3 = ' -f '
    if self.input_file_type == 'fq': cmd3 = ' -q '
    #cmd = "bowtie2 --end-to-end -a "+cmd3+" -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.08 --no-unal --omit-sec-seq --reorder"+' -p '+str(self.threads)+' -x '+self.index_file+' -U '+self.input_file
    cmd = "bowtie2 --end-to-end -a "+cmd3+" -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.12 --no-unal --reorder"+' -p '+str(self.threads)+' -x '+self.index_file+' -U '+self.input_file
    cmd2 = ' | '+self.samtools_path+' view -Sb -'
    if self.output_type == 'bam':
      cmd += cmd2
    #sys.stderr.write("\n"+cmd+"\n")
    subprocess.call(cmd,shell=True,stdout=self.output_handle,stderr=self.error_output_handle)
    return

  def set_samtools_path(self,samtools_path):
    if samtools_path == 'samtools': 
      self.samtools_path = samtools_path
    elif os.path.isfile(samtools_path): self.samtools_path = samtools_path
    else:
      sys.stderr.write("ERROR: no file in path "+samtools_path+"\n")
      sys.exit()
    return

  def set_output_type(self,output_type):
    if output_type != 'sam' and output_type != 'bam':
      sys.stderr.write("ERROR: unsupported output type "+output_type+"\n")
      sys.exit()
    self.output_type = output_type
    return

  def set_input_type(self,input_file_type):
    if input_file_type != 'fa' and input_file_type != 'fq':
      sys.stderr.write("ERROR: unsupported output type "+output_type+"\n")
      sys.exit()
    self.input_file_type = input_file_type
    return

class GenericAlignerIndexBuilder:
  def __init__(self,aligner,input_file):
    self.aligner_choice = aligner
    self.input_file = input_file
    self.input_file_type = 'fa'
    self.error_output_handle = open('/dev/null','w')
    self.threads = cpu_count()
    return

  def execute(self,output_name):
    self.output_name = output_name
    if self.aligner_choice == 'hisat':
      self.execute_hisat_build_index()
    elif self.aligner_choice == 'bowtie2':
      self.execute_bowtie2_build_index()
    return

  def execute_hisat_build_index(self):
    cmd = 'hisat-build '+self.input_file+' '+self.output_name
    subprocess.call(cmd,shell=True,stdout=self.error_output_handle,stderr=self.error_output_handle)
    return

  def execute_bowtie2_build_index(self):
    cmd = 'bowtie2-build '+self.input_file+' '+self.output_name
    subprocess.call(cmd,shell=True,stdout=self.error_output_handle,stderr=self.error_output_handle)
    return

  def set_aligner_choice(self,aligner):
    self.aligner_choice = aligner
    return
