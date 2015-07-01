#!/usr/bin/python
import argparse, subprocess, sys,os, re

# Pre: sequences separated by newlines input either by file or STDIN
#      optionally, tempdir can be specified for where temp files will go
#      (default is /tmp), and --sort_mem can be specified as the maximum
#      memory to use.
# Post: fasta file sorted by sequence with each sequence unique and each
#       name also unique output to STDOUT or file as defined.
# 
# Modifies:  This script calls itself in a weird way.  If it gives itself the
#            --level 1 flag, it will parse the STDIN as the output from uniq -c
#            and put the output fasta in the appropriate output.

def main():
  parser = argparse.ArgumentParser(description="Accept an input of sequences for a file or STDIN and output a fasta of unqiely named unique sequences")
  parser.add_argument('-i','--input',default='-',help='INPUTFILE or "-" for default of STDIN')
  parser.add_argument('-o','--output',default='-',help='OUTPUTFILE or "-" for default of STDOUT')
  parser.add_argument('--tempdir',default='/tmp',help='FOLDERNAME path to write temporary files')
  parser.add_argument('--sort_mem',type=int,help="INT amount of memory to use in sort")
  parser.add_argument('--level',type=int,default=0)
  args = parser.parse_args()

  bindir = os.path.dirname(os.path.realpath(__file__)).rstrip('/')
  inf = sys.stdin
  of = sys.stdout
  if args.output != '-': of = open(args.output,'w')
  if args.input != '-': inf = open(args.input)
  if args.level==0: # mode thats run when the script is intially run by the user
    cmd1 = "sort -T "+args.tempdir
    if args.sort_mem:
      cmd1 += " -S "+str(args.sort_mem)
    cmd1 += " | uniq -c | "+bindir+"/seq2uniqfasta.py --level 1 --output "+args.output
    p1 = subprocess.Popen(cmd1,stdin=subprocess.PIPE,shell=True)
    for line in inf:
      p1.stdin.write(line)
    p1.communicate()
  elif args.level==1: #mode thats run when the script calls itself
    z = 0
    for line in sys.stdin:
      z+=1
      line = line.rstrip()
      m = re.match('^\s*(\d+)\s+(\S+)$',line)
      if not m: 
        sys.stderr.write("ERROR unable to process uniq -c line "+line)
      of.write(">"+str(z)+"_"+m.group(1)+"\n"+m.group(2)+"\n")
  of.close()
        
if __name__=="__main__":
  main()
