# -*- coding: utf-8 -*-
"""
Created on Mon May 26 11:23:27 2014

@author: tim
"""

from maybrain import allen
from os import path
from glob import glob

probeNumbers = ["1014908", "1010837", "1012542"]
metricDict = {"degree":"brain_d2_degree_wt.txt",
              "pcCentNm":"brain_d2_pcCentNmWt_wt.txt",
              "le":"brain_local_0.03_d2_le.txt"}
#metricDict = {"degree":"brain_d2_degree_wt.txt"}
              
nodesToExclude = [28, 303, 355, 339, 131, 250, 491,205, 423, 140, 434, 142,
                  235, 244, 493, 21, 26, 232, 76, 234, 422]

diags = ["Control"] #"Controltest", "Controlretest"] # "Control", 

for d in diags:
    a = allen.multiSubj(path.join(d, "mean_2d_500.txt"), nodesToExclude=nodesToExclude)
    a.comparison()
    a.writeXMatrix(outFile=path.join(d,"Xmatrix.csv"))
    
    subjList = [v for v in glob(path.join(d, "*")) if path.isdir(v) ]
    a.writeYMatrixIndividuals(metricDict=metricDict, subjList=subjList, outFile=path.join(d,"YmatrixInd.csv"))
    a.writeYMatrixGroup(metricDict=metricDict, outFile=path.join(d,"YmatrixGroup.csv"))
    
    f = open(path.join(d,"Xmatrix.csv"), "r")
    lines = [v for v in f.readlines()]
    f.close()
        
    # remove probes for unknown genes
    f = open(path.join(d, "Xmatrix.csv"), "w")
    for l in lines:
        if l[:2]=="A_" or l[:4]=="CUST" or l[:2]=="na":
            f.writelines("# "+l)
        else:
            f.writelines(l)
    f.close()
