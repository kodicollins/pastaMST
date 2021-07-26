'''
Author: Kodi Taraszka
Email: koditaraszka@ucla.edu
This script handles the parser and i/o
'''

import argparse

class InOut():
  
  def __init__(self):
    self.indir = ''
    self.outdir = ''
    self.opal = ''
    self.final = ''
    self.MAX = False
    self.method = ''
    self.consec = False
    self.split = ''
    self.step = ''
    #print('load')
    #self.define_parser()

  
  # Add command line arguments
  def define_parser(self):
    print('hello')
    parser = argparse.ArgumentParser(description = "This program generates a minimum (or maximum) spanning tree of PASTA sub-alignment")
  
    parsargs = parser.add_argument_group("Arguments")
    parsargs.add_argument('-i', '--input', dest = 'indir', required = False, default = 'pastajob',
      help = "Full path to directory of temporary subdirectories of PASTA (pastajob). Default: pastajob")
    parsargs.add_argument('-o', '--output', dest = 'outdir', required = False, default = 'opal',
      help = "Full path to directory where OPAL alignments are stored. Default: opal")
    parsargs.add_argument('-j', '--opal', dest = 'opal', required = False, default = 'opal.jar',
      help = "Full path to OPAL jar. Default: opal.jar")
    parsargs.add_argument('-f', '--final', dest = 'final', required = False, default = 'alignedMST.fasta',
      help = "Full path to output final alignment of PASTA after using MST of sub-alignments. Default = alignedMST.fasta")
    parsargs.add_argument('-m', '--max', dest = 'MAX', action = 'store_true', default = False, required = False,
      help = "Flag for making the algorithm a maximum spanning tree instead of minimum spanning tree. Default: False")
    parsargs.add_argument('-x', '--method', dest = 'method', required = False, default = 'median',
      help = "Method for scoring alignment. Options = {median, mean, max, min}, Default: median")
    parsargs.add_argument('-c', '--consec', dest = 'consec', required = False, action = 'store_true', default = False,
      help = "Flag for changing from overall number of indels to longest consecutive indel. Default: False")
    parsargs.add_argument('-t', '--step', dest = 'step', required = False, default = '2',
      help = "PASTA step to align using MST. Default: 2 (min: 0, max: 2 if using default PASTA)")
    parsargs.add_argument('-s', '--split', dest = 'split', required = False, default = 'centroid',
      help = "PASTA splitting method, Options: {centroid, longest}, Default: centroid") 
    # Could be done differently but won't
    
    args = parser.parse_args()
    self.indir = args.indir
    if not self.indir.endswith('/'):
      self.indir += '/'
    self.outdir = args.outdir
    if not self.outdir.endswith('/'):
      self.outdir += '/'
    self.opal = args.opal
    self.final = args.final
    self.MAX = args.MAX
    self.method = args.method
    self.consec = args.consec
    self.step = args.step
    self.split = args.split
