# -*- coding: utf-8 -*-
"""
Created on Mon May 26 11:23:27 2014

@author: tim
"""

from maybrain import allen
from os import path
from glob import glob

#input association matrix
inFile = "Control/mean_500.txt"

# output matrix files
outXMat = "Xmatrix.csv"
outYMatInd = "YmatrixInd.csv"
outYmatGp = "YmatrixGroup.csv"

# Spatial information for the association matrix, with 4 columns
# 1st column contains nodes names in the same order as the association matrix
# columns 2,3,4 contain x,y,z coordinates respectively, separated by any white space delimiter
spatialFile = "parcel_500.txt"

# get a list of subjects and set a path to where the group mean data is contained
subjList = [v for v in glob(path.join("Control", "*")) if path.isdir(v) ] # list of directories containing individual subject graph metrics
pathGroup = "Control"  # path to the directory containing group level graph metrics (ie from a mean association matrix)

# probes of interest from the Allen data, if this is blank then all probes are included
probeNumbers = ["1014908"] #, "1010837", "1012542", "1010655", "1020182"]
#probeNumbers = []

# Metrics of interest with associated filename
metricDict = {"degree":"brain_d2_degree_wt.txt",
              "pcCentNm":"brain_d2_pcCentNmWt_wt.txt",
              "le":"brain_local_3_d2_le.txt"}

#nodesToExclude = [28, 303, 355, 339, 131, 250, 491,205, 423, 140, 434, 142,
#                  235, 244, 493, 21, 26, 232, 76, 234, 422]
nodesToExclude=[]

# set up the maybrain object with the group mean association matrix
a = allen.multiSubj(inFile, nodesToExclude=nodesToExclude, delim=",",
                    symmetrise=True, convertMNI=True, mirror=True)

# compare nodes to get the spatially neareest nodes from the Allen data and association matrix
a.comparison()

# write the X matrix
a.writeXMatrix(outFile=outXMat, probeNumbers=probeNumbers, tempMatName="tempMat.txt", sd=True)
#
## write out two Y matrices, one including individual subject data, the other for hte group level data
#a.writeYMatrixIndividuals(metricDict=metricDict, subjList=subjList, outFile=outYMatInd)
#a.writeYMatrixGroup(metricDict=metricDict, outFile=outYmatGp)
#
## remove probes for unknown genes from the X matrix file
#f = open(outXMat, "r")
#lines = [v for v in f.readlines()]
#f.close()
#    
#f = open(outXMat, "w")
#for l in lines:
#    if l[:2]=="A_" or l[:4]=="CUST" or l[:2]=="na":
#        f.writelines("# "+l)
#    else:
#        f.writelines(l)
#f.close()
