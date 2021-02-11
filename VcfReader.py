# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:52:50 2015

@author: hilmar
"""

import bz2file, gzip, fileinput, sys
import pysam 

class VcfReader:
    """
    Read comfortably from VCF style files with main focus on chr, start, ref and alt fields.
    Note that this API uses 1-based coordinates with both start and end included in the interval. 
    PySam API uses 0-based half-open intervals, so we have to convert internally.
    """

    def __init__(self, input_file):
        self.filename = input_file
        self.indexed = False
        
        if input_file.strip()=="-":
            ifile=sys.stdin
        elif input_file.endswith(".bz2"):
            try:        
                ifile = bz2file.BZ2File(input_file,"r", buffering=0)  
            except Exception, e:
                raise e
        elif input_file.endswith(".gz") or input_file.endswith(".bgz"):
            # try to open the file with Tabix
            try:
                ifile = pysam.Tabixfile(input_file, parser=pysam.asVCF())
                self.indexed = True
            except Exception, e:
                try:
                    ifile=gzip.GzipFile(input_file,"r")
                except Exception, e:
                    raise e
        else:
            try:
                ifile=open(input_file,"r", buffering=0)
            except Exception, e:
                raise e

        self.file = ifile
        self.last_seq = None
        self.last_start = None
        self.last_end = None
        self.__cache = None
        
    def get_block(self, seq=None, start=-1, end=-1, blocksize_max=1e6):
        """
        Retrieves a block from file and returns it as a generator
        """
        if self.indexed:
            if seq is None: # read complete file sequence by sequence, i.e. first block is first seq
                def result():
                    for s in self.file.contigs:
                        try: 
                            r = self.file.fetch(s)
                        except KeyError:
                            r = []
                        yield r
                return result()
            else:
                if not seq in self.file.contigs:
                    return []
                try: 
                    if start>0 and end>0:
                        r = self.file.fetch(seq, start-1, end)
                    else:
                        r = self.file.fetch(seq)
                except KeyError:
                    r = []
                except Exception, e:
                    print "Error while trying to accessing VCF: %s" % e.message
                    r = []
              
                return r
        else:
            raise Exception,"Non-indexed access currently not supported. Use Tabix indexed input files."
    
    def get_contigs(self):
        if self.indexed:
            return self.file.contigs
        else:
            return None

if __name__ == '__main__':
    v = VcfReader("./SRR1171905_filtered.vcf.gz")
    x = [e for e in v.get_block("19")]
    print len(x)
    for c in v.get_block():
        x = [e for e in c]
        print x[0].contig, len(x)
        