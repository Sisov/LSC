#!/usr/bin/python

import sys
import subprocess
import datetime
import random
from SequenceBasics import GenericFastaFileReader

################################################################################
# Debug flags
printSCD = False

################################################################################

# This class provides a means to convert a nav file into a map file.
# The formats of these files are not yet defined in the literature here.
# Nevertheless, this class will let you read a nav file and output results
# entry by entry
# Pre: A sorted nav file (sorted by long read name, field 2)
#      An LR cps fasta filename
#      An LR idx filename
#      LR readnames filename
#      short read coverage threshold (integer)
# Usage:
#   Create a new NavToMapFactory object
#   call read_entry until its empty
#   close it.
# Post: read_entry() outputs the line

class NavToMapFactory:
  def __init__(self,nav_filename,LR_cps_filename,LR_idx_filename,LR_readnames, short_read_coverage_threshold):  
    self.SR_cvrg_threshold = short_read_coverage_threshold
    # Read the long reads cps file
    gfr = GenericFastaFileReader(LR_cps_filename)
    self.LR_cps_dict={}
    while True:
      entry = gfr.read_entry()
      if not entry: break
      self.LR_cps_dict[entry['name']]=entry['seq']

    # Read the long reads idx file
    self.LR_idx_dict={}
    with open(LR_idx_filename,'r') as LR_idx_file:
      for line in LR_idx_file:
        fields=line.strip().split('\t')
        self.LR_idx_dict[fields[0]] = '\t'.join(fields[1:])

    # Read the long read names
    self.LR_name_dict={}
    with open(LR_readnames) as LR_name_file:
      for line in LR_name_file:
        fields=line.strip().split('\t')
        self.LR_name_dict[fields[0]] = fields[1]

    #Initialize the begining
    with open(nav_filename ,'r') as nav_cvrg_file:   # This used to compute coverage
      line_list = nav_cvrg_file.readline().strip().split('\t')

    self.LR_name = line_list[1]
    self.LR_coverage_list  = [0] * len(self.LR_cps_dict[self.LR_name])
    self.LR_uniq_coverage_list = [0] * len(self.LR_cps_dict[self.LR_name])
    self.SR_rpt = int(line_list[0].split('_')[1])
    self.pos = int(line_list[2])
    self.SR_len = len(line_list[3])
    self.LR_coverage_list[self.pos:(self.pos+self.SR_len)] = [(i + self.SR_rpt) for i in self.LR_coverage_list[self.pos:(self.pos+self.SR_len)]]
    self.LR_uniq_coverage_list[self.pos:(self.pos+self.SR_len)] = [(i + 1) for i in self.LR_uniq_coverage_list[self.pos:(self.pos+self.SR_len)]]
    self.nav_cvrg_file = open(nav_filename ,'r')
    self.prev_LR_name = self.LR_name
    self.enable_printSCD = False
    self.buffer = []

  def read_entry(self):
    while True: # Stay in until we get an entry or an end of file
      line = self.nav_cvrg_file.readline()
      if not line and len(self.buffer) == 0:
        return None
      elif not line: # End of file but flush the buffer
        oline = write_2_LR_SR_map_file(self.buffer, self.LR_coverage_list, self.LR_uniq_coverage_list,self.SR_cvrg_threshold,self.LR_cps_dict,self.LR_idx_dict)
        self.buffer = []
        return oline
      # Not the end of the file
      line_list=line.strip().split('\t')

      self.LR_name = line_list[1]
      oline = None
      if (self.prev_LR_name != self.LR_name):
            oline = write_2_LR_SR_map_file(self.buffer, self.LR_coverage_list, self.LR_uniq_coverage_list,self.SR_cvrg_threshold,self.LR_cps_dict,self.LR_idx_dict)
            #LR_SR_mapping_file.write(oline+"\n")
            
            self.LR_coverage_list  = [0] * len(self.LR_cps_dict[self.LR_name])
            self.LR_uniq_coverage_list  = [0] * len(self.LR_cps_dict[self.LR_name])
            self.prev_LR_name = self.LR_name
            self.buffer = []
      self.buffer.append(line)
      self.SR_rpt = int(line_list[0].split('_')[1])
      self.pos = int(line_list[2])
      self.SR_len = len(line_list[3])
      self.LR_coverage_list[self.pos:(self.pos+self.SR_len)] = [(i + self.SR_rpt) for i in self.LR_coverage_list[self.pos:(self.pos+self.SR_len)]]
      self.LR_uniq_coverage_list[self.pos:(self.pos+self.SR_len)] = [(i + 1) for i in self.LR_uniq_coverage_list[self.pos:(self.pos+self.SR_len)]]
      if oline: return oline

  def enable_printSCD(self,temp_foldername):
    self.enable_printSCD = True
    temp_foldername = temp_foldername.rstrip('/')+'/'
    self.LR_SR_coverage_file = open(temp_foldername + "LR_SR.scd",'w')
    self.LR_uSR_coverage_file = open(temp_foldername + "LR_SR.uscd",'w')
    self.LR_SR_coverage_selected_file = open(temp_foldername + "LR_SR.scd.selected",'w')
    self.LR_uSR_coverage_selected_file = open(temp_foldername + "LR_SR.uscd.selected",'w')

  def close(self):
    if self.enable_printSCD:
      self.LR_SR_coverage_file.close()
      self.LR_uSR_coverage_file.close()
      self.LR_SR_coverage_selected_file.close()
      self.LR_uSR_coverage_selected_file.close()
    
# This function does an unknown conversion of a nav to a map file.
def write_2_LR_SR_map_file(buffer, LR_coverage_list, LR_uniq_coverage_list,SR_cvrg_threshold, LR_cps_dict,LR_idx_dict):
    
    LR_name = buffer[0].rstrip().split("\t")[1]
    # Store LR-SR SCD
    if (printSCD):
        LR_SR_coverage_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_coverage_list]
        LR_SR_coverage_file.write(",".join(LR_coverage_str_list) + "\n")
        
        LR_uSR_coverage_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_uniq_coverage_list]
        LR_uSR_coverage_file.write(",".join(LR_coverage_str_list) + "\n")
        
        LR_coverage_str_list_temp = [0] * len(LR_coverage_str_list)
        LR_uniq_coverage_str_list_temp = [0] * len(LR_coverage_str_list)
    
    LR_SR_list = []
    for line in buffer:
        if line[0]=='#': continue
        line_list=line.strip().split('\t')
        
        SR_rpt = int(line_list[0].split('_')[1])
        pos = int(line_list[2])
        SR_len = len(line_list[3])
        region_cvrg = min(LR_coverage_list[pos:(pos+SR_len)])
        if (SR_cvrg_threshold < 0):    # Disable the feature for negative threshold
            SR_prob = 1.1
        else:
            SR_prob = 1./region_cvrg * SR_cvrg_threshold   # Note: SR_prob could be greater than 1
        
        use_SR = False
        for SR_rpt_idx in range(SR_rpt):
            if (random.random() <= SR_prob):
                use_SR = True
                break
        if not use_SR: continue
        if (printSCD):
                LR_coverage_str_list_temp[pos:(pos+SR_len)] = [i + SR_rpt for i in  LR_coverage_str_list_temp[pos:(pos+SR_len)]]
                LR_uniq_coverage_str_list_temp[pos:(pos+SR_len)] = [i + 1 for i in  LR_uniq_coverage_str_list_temp[pos:(pos+SR_len)]]
        if (line_list[3] == "*"):
                line_list[3] = ""
        line_list_temp = [line_list[0],line_list[2]] + line_list[3:]
        LR_SR_list.append(line_list_temp)
    if (printSCD):
        LR_SR_coverage_selected_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_coverage_str_list_temp]
        LR_SR_coverage_selected_file.write(",".join(LR_coverage_str_list) + "\n")
        
        LR_uSR_coverage_selected_file.write(">" + LR_name_dict[LR_name] + "\n")
        LR_coverage_str_list = [str(i) for i in LR_uniq_coverage_str_list_temp]
        LR_uSR_coverage_selected_file.write(",".join(LR_coverage_str_list) + "\n")

    temp_SR_ls = []
    ls_SR_seq = []
    ls_SR_idx_seq = []
    for SR in LR_SR_list:
        temp_SR_ls.append(SR[0]+','+ SR[1]+','+SR[2])
        ls_SR_seq.append(SR[3])
        if (len(SR) > 4):
            ls_SR_idx_seq.append('\t'.join([SR[4], SR[5]]))
        else:
            ls_SR_idx_seq.append('\t'.join(["", ""]))   # no compression point

    input_ls= [LR_cps_dict[LR_name], LR_idx_dict[LR_name], ';'.join(temp_SR_ls),LR_name,'kinfai'.join(ls_SR_seq), 'kinfai'.join(ls_SR_idx_seq)]
    return 'yue'.join(input_ls)

if __name__ == "__main__":
  main()
