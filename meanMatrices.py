# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 23:24:39 2012

@author: tim
"""

import os,glob,csv
from os import path

fNames = glob.glob("PD/*/wave_cor_mat_level_2d_500.txt")

for i in range(500):
    if not path.exists("dir_"+str(i)):
        os.mkdir("dir_"+str(i))
    f = open(path.join("dir_"+str(i),"allValues.txt"), "wb")
    f.close()
    
for fName in fNames:
    print 'loading data from', fName
    # open file
    f = open(fName,"rb")
    reader = csv.reader(f, delimiter=" ")
  
    for i in range(500):
        g = open(path.join("dir_"+str(i),"allValues.txt"), "ab")
        writer = csv.writer(g, delimiter=",")
	l = reader.next()
        line = [v if v!="nan" else 0 for v in l]
        writer.writerow(line)
        del(line)
        g.close()
    del(reader)
    f.close()
