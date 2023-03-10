"""
Extract the resolved reads from the raw reads 

"""
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import time
from gzip import open as gzopen
import gzip
from functools import partial
import sys

if __name__ == '__main__':
      ''' 
      resolved_file = 'TEST_NV_NW_resolved.fq.gz'

      fq1 = 'TEST_R1.fq.gz'
      fq2 = 'TEST_R2.fq.gz'

      fq1_out = 'resolved_TEST_R1.fq.gz'
      fq2_out = 'resolved_TEST_R2.fq.gz'
      '''

      resolved_file = sys.argv[1]

      fq1 = sys.argv[2]
      fq2 = sys.argv[3]

      fq1_out = 'resolved_' + fq1
      fq2_out = 'resolved_' + fq2


      _openz = partial(gzip.open, mode='wt') 

      required = { r.id.split(':cseq')[0][:-3] for r in SeqIO.parse(gzopen(resolved_file,'rt'), 'fastq')  } # required readnames 
      print ('From file:', resolved_file, ' total_reads = ', len(required) )
      
      with  _openz(fq1_out) as outF1,  _openz(fq2_out) as outF2:
        for r1, r2 in zip(SeqIO.parse(gzopen(fq1,'rt'), 'fastq'), SeqIO.parse(gzopen(fq2,'rt'), 'fastq')):
            if r1.id in required:
               SeqIO.write(r1, outF1, "fastq")
               SeqIO.write(r2, outF2, "fastq")
      
      print ('Done')

