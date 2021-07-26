'''
Author: Kodi Taraszka
Email: koditaraszka@ucla.edu
This script is the main script that inherits information from IO and Methods
'''

from inout import InOut
from methods import Methods
import glob
import multiprocessing 

class Main(InOut, Methods):
  def __init__(self):
    super().__init__()
    self.nproc = multiprocessing.cpu_count()
    self.mafft_dict = {}
    #self.mafft_score = []
    self.merged_groups = {}
    self.merged = {}
    self.parent = {}
    self.rank = {}
    self.mst = []
    self.alignment = {}
    self.opal_dict = {}
    self.MST()

  def MST(self):
    self.getMafftAlignment()
    self.opalPairwise()
    sort_keys = sorted(self.opal_dict.keys(), reverse = self.MAX)
    #For MAX = TRUE, reverse = TRUE
    #TODO: what is the key?
    for k in sort_keys:
      # Loop over pairwise opal alignment
      # TODO: how does find work? 
      for tups in self.opal_dict[k]:
        if self.find(tups[1]) != self.find(tups[2]):
          self.union(tups[1], tups[2])
          self.mst.append(tups)
    self.trans()

  def score(self, opal):
    seqs = all.Alignment()
    seqs.read_file_object(opal)
    gap = list()
    if self.consec:
      for i in seqs:
        theSeq = list(seqs[i])
        longest = 0
        current = 0
        possible = len(theSeq)
        for j in seq(0,len(theSeq)):
          if theSeq[j] == "-":
            current += 1
          else:
            if current > longest:
              longest = current    
            current = 0
          if ((possible - j) + current) < longest:
            break
        gap.append(longest)
      gap.sort()
          
    else:
      for i in seqs:
        theSeq = list(seqs[i])
        gap.append(theSeq.count("-")/float(len(theSeq)))
      gap.sort()

    if self.method == 'median':
      if len(gap)%2 == 0:
        return gap[len(gap)/2]
      else:
        return (gap[len(gap)/2] + gap[len(gap)/2 + 1])/2.0

    elif self.method == 'max':
        return max(gap)
 
    elif self.method == 'min':
        return min(gap)

    else: #mean
        return (sum(gap)/float(len(gap)))

  def opalPairwise(self):
    argslist = []
    #this k1,k2 is a pairwise look at keys to mafft alignments
    for k1, k2 in itertools.combinations(self.mafft_dict, 2):
      opalName = self.outdir + str(k1)+'_'+str(k2)+'.txt'
      argslist.append((self.mafft_dict[k1][0],self.mafft_dict[k2][0], opalName))
      results = Pool(self.nproc).map(self.runOpal, argslist)
      for i in results:
        if i[0] in self.opal_dict:
          self.opal_dict[i[0]].append((i[1],i[2],i[3]))
        else:
          self.opal_dict[i[0]] = [(i[1],i[2],i[3])]

  def runOpal(self, args):
    mafft1 = str(args[0]) 
    mafft2 = str(args[1])
    output = str(args[2])
    args = ['java', '-Xmx1000m', '-jar', self.opal, '--in', mafft1, '--in2', mafft2, 
                '--out', output, '--quiet', '--align_method', 'profile']
    subprocess.call(args)
    x = self.score(output)
    return [x, output, mafft1, mafft2]


  def getMafftAlignment(self):
    keyName = 1
    #this only works for pasta and how it stores the data, it is not generalized
    for loc in glob.glob(self.indir+'/temp*/step'+self.step+'/'+self.split+'/r*/d*/*/input.aligned'):  
      self.mafft_dict[keyName] = [loc, mafft_score]
      keyName += 1

if __name__=="__main__":
  Main()
