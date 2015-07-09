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

def main():
  if len(sys.argv) >= 2:
    temp_foldername = sys.argv[1]
    nav_filename = sys.argv[2] #must be sorted by long read
    LR_filename = sys.argv[3]
    SR_cvrg_threshold = int(sys.argv[4])
    sort_max_mem = sys.argv[5]
    LR_readnames = sys.argv[6]
    output_filename = sys.argv[7]
    
  else:
    sys.stderr.write("usage: python convertNAV.py temp_foldername LR_filename sorted_nav_filename sort_max_mem\n")
    sys.stderr.write("or ./convertNAV.py temp_foldername LR_filename nav_filename Nthread sort_max_mem\n")
    sys.exit(1)

  # Read the long reads cps file
  gfr = GenericFastaFileReader(LR_filename+'.cps')
  LR_cps_dict={}
  while True:
    entry = gfr.read_entry()
    if not entry: break
    LR_cps_dict[entry['name']]=entry['seq']

  # Read the long reads idx file
  LR_idx_dict={}
  with open(LR_filename + '.idx','r') as LR_idx_file:
    for line in LR_idx_file:
      fields=line.strip().split('\t')
      LR_idx_dict[fields[0]] = '\t'.join(fields[1:])

  # Read the long read names
  LR_name_dict={}
  with open(LR_readnames) as LR_name_file:
    for line in LR_name_file:
      fields=line.strip().split('\t')
      LR_name_dict[fields[0]] = fields[1]

  LR_SR_mapping_file = open(output_filename,'w')
  if (printSCD):
    LR_SR_coverage_file = open(temp_foldername + "LR_SR.scd",'w')
    LR_uSR_coverage_file = open(temp_foldername + "LR_SR.uscd",'w')
    LR_SR_coverage_selected_file = open(temp_foldername + "LR_SR.scd.selected",'w')
    LR_uSR_coverage_selected_file = open(temp_foldername + "LR_SR.uscd.selected",'w')

  with open(nav_filename ,'r') as nav_cvrg_file:   # This used to compute coverage
    line_list = nav_cvrg_file.readline().strip().split('\t')

  LR_name = line_list[1]
  LR_coverage_list  = [0] * len(LR_cps_dict[LR_name])
  LR_uniq_coverage_list = [0] * len(LR_cps_dict[LR_name])
  SR_rpt = int(line_list[0].split('_')[1])
  pos = int(line_list[2])
  SR_len = len(line_list[3])
  LR_coverage_list[pos:(pos+SR_len)] = [(i + SR_rpt) for i in LR_coverage_list[pos:(pos+SR_len)]]
  LR_uniq_coverage_list[pos:(pos+SR_len)] = [(i + 1) for i in LR_uniq_coverage_list[pos:(pos+SR_len)]]

  prev_LR_name = LR_name
  num_lines = 0   # pre-incremented

  buffer = []
  with open(nav_filename ,'r') as nav_cvrg_file:   # This used to compute coverage
    for line in nav_cvrg_file:
      num_lines += 1
      if line[0]=='#': continue
         
      line_list=line.strip().split('\t')

      LR_name = line_list[1]
      if (prev_LR_name != LR_name):
            oline = write_2_LR_SR_map_file(buffer, LR_coverage_list, LR_uniq_coverage_list,SR_cvrg_threshold,LR_cps_dict,LR_idx_dict)
            LR_SR_mapping_file.write(oline+"\n")
            
            LR_coverage_list  = [0] * len(LR_cps_dict[LR_name])
            LR_uniq_coverage_list  = [0] * len(LR_cps_dict[LR_name])
            num_lines = 0
            prev_LR_name = LR_name
            buffer = []

      buffer.append(line)

      SR_rpt = int(line_list[0].split('_')[1])
      pos = int(line_list[2])
      SR_len = len(line_list[3])
      LR_coverage_list[pos:(pos+SR_len)] = [(i + SR_rpt) for i in LR_coverage_list[pos:(pos+SR_len)]]
      LR_uniq_coverage_list[pos:(pos+SR_len)] = [(i + 1) for i in LR_uniq_coverage_list[pos:(pos+SR_len)]]

  num_lines += 1
  if len(buffer) > 0:
    oline = write_2_LR_SR_map_file(buffer, LR_coverage_list, LR_uniq_coverage_list, SR_cvrg_threshold,LR_cps_dict,LR_idx_dict)
    LR_SR_mapping_file.write(oline+"\n")
  LR_SR_mapping_file.close()
  if (printSCD):
    LR_SR_coverage_file.close()
    LR_uSR_coverage_file.close()
    LR_SR_coverage_selected_file.close()
    LR_uSR_coverage_selected_file.close()
    
  sys.stderr.write("Done with generating LR_SR.map file\n")

def write_2_LR_SR_map_file(buffer, LR_coverage_list, LR_uniq_coverage_list,SR_cvrg_threshold, LR_cps_dict,LR_idx_dict):
    
    LR_name = buffer[0].rstrip().split("\t")[1]
    print LR_name +"\t long read"
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
    #LR_SR_mapping_file.write('yue'.join(input_ls)+'\n')
    return 'yue'.join(input_ls)

if __name__ == "__main__":
  main()
