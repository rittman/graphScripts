# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 23:50:09 2012

@author: tim
"""

import numpy as np
from glob import glob
from os import chdir,path,getcwd

mainDir = getcwd()

for d in glob("dir*"):
  chdir(d)
  # load file
  arr = np.loadtxt("allValues.txt", delimiter=",")

  # calculate means
  meanarr = np.mean(arr, axis=0)
  np.savetxt("meanValues.txt", meanarr[None], delimiter=",")
  chdir(mainDir)
