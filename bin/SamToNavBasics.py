import sys, os, re, numpy
from SequenceBasics import GenericFastaFileReader

#Generate a novoalign native format file from a sam format file

rev_cmplmnt_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
rev_cmplmnt_bases = rev_cmplmnt_map.keys()
def reverse_complement(seq):
    
    output_seq = ''
    len_seq = len(seq)
    for i in range(len_seq):
        if (seq[len_seq - 1 - i] in rev_cmplmnt_bases):
            output_seq += rev_cmplmnt_map[seq[len_seq - 1 - i]]
        else:
            print "Err: Unexpected base in short read sequence: " + seq
            output_seq += seq[len_seq - 1 - i]
        
    return output_seq
    

def get_SR_sequence(SR_file, SR_idx_file, SR_seq_readname):
    read_name = "invalid_read_name"
    while (read_name != SR_seq_readname):
        read_name = (SR_file.readline().strip())
        if (not read_name):
            print "Err: HC SR sequence not found."
            exit(1) 
        if (read_name[0] != '>'):
            continue    # unexpected string
        read_name = read_name[1:]
        SR_seq = SR_file.readline().strip().upper()
        SR_idx_seq = SR_idx_file.readline().strip().split('\t')
        if(SR_idx_seq[0] != read_name):
            print "Err: SR.fa.cps and SR.fa.idx do not match."
            exit(1) 
        SR_idx_seq = '\t'.join(SR_idx_seq[1:])
    return [SR_seq, SR_idx_seq]

class SamToNavFactory:
  #SR_file (query) and SR_idx_file must go in the same order as 
  #the sam file.
  def __init__(self):
    self.error_rate_threshold = 100
    self.one_line_per_alignment = False # BWA will put them all on one line
    self.current_query = "invalid read name"
    return

  def set_error_rate_threshold(self,error_rate_threshold):
    self.error_rate_threshold = error_rate_threshold
    return

  def set_one_line_per_alginment(self,one_line_per_alignment):
    self.one_line_per_alignment = one_line_per_alignment
    return

  def initialize_compressed_query(self,SR_filename,SR_idx_filename):
    self.SR_file = open(SR_filename)
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
    SR_seq = ""
    SR_seq_rvs_cmplmnt = ""
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
        [SR_seq, SR_idx_seq] = get_SR_sequence(self.SR_file, self.SR_idx_file, SR_name)
        SR_seq_rvs_cmplmnt = reverse_complement(SR_seq)
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
        

    read_seq_len = len(SR_seq)
    ostrings = []
    for align_str in align_list:
        err_state = False
        fields = align_str.split(',')
        
        ref_seq = self.LR_seq[fields[0]]
        ref_seq_len = len(ref_seq)
        if (fields[1][0] == '-'):     # Check if seq is reversed complement
            read_seq = SR_seq_rvs_cmplmnt
            pseudo_SR_name = "-" + SR_name
        else:
            read_seq = SR_seq
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
                    read_subseq = numpy.array(list(read_seq[read_idx:(read_idx + subseq_len)]))
                    ref_subseq = numpy.array(list(ref_seq[ref_idx:(ref_idx + subseq_len)]))
                    mut_indices = numpy.where((ref_subseq == read_subseq) == False)[0].tolist()
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
                ostrings.append('\t'.join([pseudo_SR_name, fields[0], fields[1], ' '.join(diff_list),SR_seq, SR_idx_seq]))
    if len(ostrings) == 0: return None
    return ostrings
  def close(self):
    self.SR_file.close()
    self.SR_idx_file.close()
