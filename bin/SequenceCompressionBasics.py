class HomopolymerCompressionFactory:
  def __init__(self):
    self.MinNonN=None
    self.MaxN=None
  def set_MinNonN(self,input_MinNonN):
    self.MinNonN = input_MinNonN
  def set_MaxN(self,input_MaxN):
    self.MaxN = input_MaxN

  # Pre: Take a nucleotide sequence (can contain N's)
  # Post: Return None if violates MinNonN or MaxN criteria.
  #       Return 1. compressed sequence
  #              2. position array
  #              3. length array
  def compress(self,seq):
    seq = seq.upper() #Make sure we are in an all uppercase situation

    # Return values
    cps_seq = seq[0] # Compressed sequence
    pos_ls=[]
    len_ls = []
    n_count = 0
    
    # Loop carry values
    repeat_count=0
    ref_s = seq[0]
    i=0
    for s in seq:
      if not ref_s == s:
        if (ref_s == 'N'):
          n_count += 1
        if repeat_count>1:
          len_ls.append(repeat_count)
          pos_ls.append(i) 
        cps_seq = cps_seq + s
        repeat_count=1
        ref_s = s
        i+=1
      else:
        repeat_count+=1
    if (ref_s == 'N'):
      n_count += 1
    if repeat_count>1:
      len_ls.append(repeat_count)
      pos_ls.append(i) 

    if self.MinNonN: 
      if len(cps_seq)-n_count < self.MinNonN: return None
    if self.MaxN:
      if n_count > MaxN: return None
    hpc = {}
    return [cps_seq, ','.join([str(x) for x in pos_ls]), ','.join([str(x) for x in len_ls])]

