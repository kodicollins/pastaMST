'''
Author: Kodi Taraszka
Email: koditaraszka@ucla.edu
This script handles the actual implementation.
'''

import pasta.alignment as all

class Methods():
  def __init__(self):
      self.merge_groups = {}
      self.merged = {}
      self.parent = {}
      self.rank = {}
      self.mst = []

  # NO IDEA what/why for makeSet, find, union
  def makeSet(self, v):
    self.parent[v] = v
    self.rank[v] = 0

  def find(self, v):
    if v not in self.parent:
      self.makeSet(v)
    elif self.parent[v] != v:
      self.parent[v] = self.find(parent[v])
    return self.parent[v]

  def union(self, v1, v2):
    root1 = self.find(v1)
    root2 = self.find(v2)
    if root1 != root2:
      if self.rank[root1] > self.rank[root2]:
        self.parent[root2] = root1
      elif self.rank[root2] < self.rank[root1]:
        self.parent[root1] = root2
      else 
          self.parent[root1] == root2
          self.rank[root2] += 1

  def findHome(self, temp)
    for i in self.merged_groups:
      if temp[1] in self.merged_groups[i] or temp[2] in self.merged_groups[i]:        
        #If it's in one, add it to both
        if temp[1] not in self.merged_groups[i]:
          self.merged_groups[i].append(temp[1])
        if temp[2] not in self.merged_groups[i]:
          self.merged_groups[i].append(temp[2])
        self.seq.read_file_object(temp[0])
        self.alignment[i].merge_in(self.seq)
        self.merged[i].append(temp[0])
        return i
      return -1 

  def round2(self):
    for k1, k2 in itertools.combinations(merged_groups, 2):
      for i in merged_groups[k1]:
        for j in merged_groups[k2]:
          if i == j:
              self.merged_groups[k2].remove(j)
              self.merge_groups[k1].extend(self.merged_groups[k2])
              self.merged[k1].extend(self.merged[k2])
              self.alignment[k1].merge_in(self.alignment[k2])
              self.merged.pop(k2)
              self.merged_groups.pop(k2)
              self.alignment.pop(k2)
              return
    return

  def trans(self)
    count = 0
    while len(self.mst) > 0:
      temp = self.mst[0]
      self.mst.pop(0)
      bestOpal = self.findHome(temp)
      if bestOpal == -1:
        self.merged[count] = [temp[0]]
        self.merged_groups[count] = [temp[1]]
        self.merged_groups[count] = [temp[2]]
        self.seq.read_file_object(temp[0])
        self.alignment[count] = self.seq
        #does count always iterate or just here, here for now
        count += 1
      self.round2()
    return

