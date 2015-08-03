import sys, re
from SequenceBasics import GenericFastaFileReader, rc

#Generate a novoalign native format file from a sam format file

### get_SR_sequence ######
#
# Read the short read sequence (compressed version) from a file that is ordered
# by short read name
#
# Pre: SR_file - an open generic fasta file reader to the homopolymer compressed fasta file
#      SR_idx_file - an open file handle to the idx file from the homopolymer compressed sequence
#      SR_seq_readname - the name of the sequence to retrieve
# Post:  Return the nucleotide sequence (compressed) and the compression scheme (2nd field onward from idx file) in a tsv format.
# Modifies: Will step forward in the file handles and fasta file reader
#           Therefore the file must be ordered so that each requested read comes sequentially later in the file
def get_SR_sequence(SR_file, SR_idx_file, SR_seq_readname):
  read_name = "invalid_read_name"
  SR_seq = None
  SR_idx_seq = None
  while (read_name != SR_seq_readname):
    entry = SR_file.read_entry()
    if not entry:
      sys.stderr.write("Err: HC SR sequence not found.\n")
      sys.exit() 
    read_name = entry['name']
    SR_seq = entry['seq'].upper()
    SR_idx_seq = SR_idx_file.readline().strip().split('\t')
    if SR_idx_seq[0] != read_name:
      sys.stderr.write("Err: SR.fa.cps and SR.fa.idx do not match.\n")
      sys.exit()
    SR_idx_seq = '\t'.join(SR_idx_seq[1:])
  return [SR_seq, SR_idx_seq]

class SamToNavFactory:
  #SR_file (query) and SR_idx_file must go in the same order as 
  #the sam file.
  def __init__(self):
    self.error_rate_threshold = 100
    self.one_line_per_alignment = False # BWA will put them all on one line
    self.current_query = "invalid read name"
    self.SR_seq = None
    self.SR_idx_seq = None
    self.SR_seq_rvs_cmplmnt = None
    self.line_count = 0
    self.thread_count = 1
    self.thread_index = 1
    return
  def set_thread_count(self,thread_count):
    self.thread_count = thread_count
  def set_thread_index(self,thread_index):
    self.thread_index = thread_index
  def set_error_rate_threshold(self,error_rate_threshold):
    self.error_rate_threshold = error_rate_threshold
    return

  def set_one_line_per_alginment(self,one_line_per_alignment):
    self.one_line_per_alignment = one_line_per_alignment
    return

  def initialize_compressed_query(self,SR_filename,SR_idx_filename):
    self.SR_file = GenericFastaFileReader(SR_filename)
    self.SR_idx_file = open(SR_idx_filename)
    return

  def initialize_target(self,LR_filename):
    gfr = GenericFastaFileReader(LR_filename)
    self.LR_seq = {}
    while True:
      entry = gfr.read_entry()
      if not entry: break
      self.LR_seq[entry['name']] = entry['seq']
    gfr.close()
    return self.LR_seq

  def sam_to_nav(self,line):
    self.line_count += 1
    if self.line_count % self.thread_count != self.thread_index-1:
      return None
    if not self.LR_seq:
      sys.stderr.write("ERROR initialize target first\n")
    if not self.SR_file:
      sys.stderr.write("ERROR initialize query first\n")
    
    if (line[0] == '@'):
        return None
    
    line_fields = line.strip().split('\t')
    cigar = line_fields[5]
    if ((cigar == '*') or (cigar == '.')):
        return None
    
    SR_name = line_fields[0]
    if (SR_name != self.current_query):
        [self.SR_seq, self.SR_idx_seq] = get_SR_sequence(self.SR_file, self.SR_idx_file, SR_name)
        self.SR_seq_rvs_cmplmnt = rc(self.SR_seq)
        self.current_query = SR_name
    if (int(line_fields[1]) & 0x10):     # Check if seq is reversed complement
        line_fields[3] = '-' + line_fields[3]
    else:
        line_fields[3] = '+' + line_fields[3]
    
    align_list = [','.join([line_fields[2], line_fields[3], line_fields[5], str(0)])]
    
    if (not self.one_line_per_alignment):   # BWA reports all alignment per read in one line
        multi_align_str = ','.join([line_fields[2], line_fields[3], line_fields[5], str(0)]) + ';'
        for fields_idx in range(11, len(line_fields)):
            if (line_fields[fields_idx][0:5] == 'XA:Z:'):
                multi_align_str += line_fields[fields_idx][5:]
                break
        align_list =  multi_align_str[:-1].split(';')
        

    read_seq_len = len(self.SR_seq)
    ostrings = []
    for align_str in align_list:
        err_state = False
        fields = align_str.split(',')
        
        ref_seq = self.LR_seq[fields[0]]
        ref_seq_len = len(ref_seq)
        if (fields[1][0] == '-'):     # Check if seq is reversed complement
            read_seq = self.SR_seq_rvs_cmplmnt
            pseudo_SR_name = "-" + SR_name
        else:
            read_seq = self.SR_seq
            pseudo_SR_name = SR_name
        fields[1] = fields[1][1:]
        read_idx = 0
        sub_ref_idx =  1  # 1-offset address
        ref_idx = int(fields[1]) - 1   # convert to 0-offset address
        diff_list = []
        cigar_list = re.split('(M|I|D)', fields[2])
        num_err = 0
        for idx in range(1, len(cigar_list), 2):
            if (cigar_list[idx - 1].isdigit()):
                if (cigar_list[idx] == 'M'):
                    subseq_len = int(cigar_list[idx - 1])
                    if ((read_idx + subseq_len > read_seq_len) or
                         (ref_idx + subseq_len > ref_seq_len)):
                        err_state = True
                        break
                    read_subseq = list(read_seq[read_idx:(read_idx + subseq_len)])
                    ref_subseq = list(ref_seq[ref_idx:(ref_idx + subseq_len)])
                    mut_indices = [x for x in range(len(read_subseq)) if read_subseq[x] != ref_subseq[x]]
                    for mut_idx in mut_indices:
                        if (read_subseq[mut_idx] != "N"):
                            diff_list += [str(sub_ref_idx + mut_idx) + ref_subseq[mut_idx] + '>' + read_subseq[mut_idx]]
                    read_idx += subseq_len
                    ref_idx += subseq_len
                    sub_ref_idx += subseq_len
                    num_err += len(mut_indices)
                elif (cigar_list[idx] == 'I'):
                    subseq_len = int(cigar_list[idx - 1])
                    if (read_idx + subseq_len > read_seq_len):
                        err_state = True
                        break
                    insert_str = re.sub(r'N|n', '', read_seq[read_idx:(read_idx + int(cigar_list[idx - 1]))])
                    if (insert_str != ""):
                        diff_list += [str(sub_ref_idx) + '+' + read_seq[read_idx:(read_idx + subseq_len)]]
                    read_idx += subseq_len
                    num_err += subseq_len
                elif (cigar_list[idx] == 'D'):
                    subseq_len = int(cigar_list[idx - 1])
                    if (ref_idx + subseq_len > ref_seq_len):
                        err_state = True
                        break
                    for del_idx in range(subseq_len):
                        diff_list += [str(sub_ref_idx + del_idx) + '-' + ref_seq[ref_idx + del_idx]]
                    ref_idx += subseq_len
                    sub_ref_idx += subseq_len
                    num_err += subseq_len
            else:
                #print 'Err in cigar: tag : ' + line
                err_state = True
                break
        if ((cigar_list[-1] == '') and (not err_state)):
            if (len(diff_list) == 0):
                diff_list.append("*")
            err_rate = (100 * num_err) / read_seq_len
            if (err_rate <= self.error_rate_threshold):
                ostrings.append('\t'.join([pseudo_SR_name, fields[0], fields[1], ' '.join(diff_list),self.SR_seq, self.SR_idx_seq]))
    if len(ostrings) == 0: return None
    return ostrings
  def close(self):
    self.SR_file.close()
    self.SR_idx_file.close()
