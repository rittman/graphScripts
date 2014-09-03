# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 00:18:58 2012

@author: tim
"""

from maybrain import mayBrainTools as mbt # main maybrain package
from maybrain import mayBrainExtraFns as extras # functions to write results and other bits
import bct # used to integrate pythonic brain connectivity toolbox
import numpy as np
import community
#from cocotools import infomap
from os import remove
from metrics import metrics

edgePCCons = [v for v in range(1,11)] # threholds 1 to 10%
excludedNodes = [28, 303, 355, 339, 131, 250, 491, 205, 423, 140, 434, 142, 235,
                 244, 493, 21, 26, 232, 76, 234, 422] # nodes with insufficient coverage

dVal = "2"  # wavelet level (only for use in filename below)
adjMatFile = "wave_cor_mat_level_"+dVal+"d_500.txt" # input association matrix file
                 244, 493, 21, 26, 232, 76, 234, 422]
nP = "500"
parcelFile = "parcel_"+nP+".txt"  spatial information in 4 columns: node, x, y, z
thresholdtype = "local" # applies local thresholding using MST

delim=" "

# create brain object
a = mbt.brainObj()
appVal = False
appValT = False

# unweighted measures
# iterate through thresholds
for e in edgePCCons:
    ofb = '_'.join(["brain", thresholdtype, str(e), "d"+dVal+"_"])
    propDict = {"edgePC":str(e)} # added properties for results files

    a.importAdjFile(adjMatFile, delimiter=delim, excludedNodes=excludedNodes)
    a.localThresholding(edgePC=e)  # apply a threshold
    a.removeUnconnectedNodes()     # remove unconnected nodes
    
    degs = a.G.degree(weight='weight')
    extras.writeResults(degs, "degreeWt", ofb, append=appVal)

    a.binarise()  # binarise the graph
    a.importSpatialInfo(parcelFile)  # read spatial information
    a.weightToDistance() # convert weights to distance (for closeness centrality function)
    a.makebctmat() # create an array to be used by the brain connectivity toolbox functions
    
    ofbT = '_'.join(["brain", thresholdtype, str(e), "d"+dVal+"_"])  
   
    #### small worldness metrics ####
    degs = mbt.nx.degree(a.G)  # measure degree
    extras.writeResults(degs, "degree", ofb, propDict=propDict, append=appVal)  # write the results to a file
        
    clustCoeff = mbt.nx.average_clustering(a.G)
    extras.writeResults(clustCoeff, "clusterCoeff", ofb, propDict=propDict, append=appVal)
    del(clustCoeff)
    
    pl = mbt.nx.average_shortest_path_length(a.G)
    extras.writeResults(pl, "pl", ofb, propDict=propDict, append=appVal)
    del(pl)
    
    ge = extras.globalefficiency(a.G)
    extras.writeResults(ge, "ge", ofb, propDict=propDict, append=appVal)
    del(ge)
    
    le = extras.localefficiency(a.G)
    extras.writeResults(le, "le", ofb, propDict=propDict, append=appVal)
    del(le)

    # hub metrics
    betCent = mbt.nx.centrality.betweenness_centrality(a.G)
    extras.writeResults(betCent, "betCent", ofb, propDict=propDict, append=appVal)
    
    closeCent = mbt.nx.centrality.closeness_centrality(a.G)
    extras.writeResults(closeCent, "closeCent", ofb, propDict=propDict, append=appVal)
     
#    hs = extras.hubscore(a.G, bc=betCent, cc=closeCent, degs=degs, weighted=False)
#    extras.writeResults(hs, "hs", ofb, propDict=propDict, append=appVal)
#    del(hs, betCent, closeCent, degs)
     
    try:
        eigCent = mbt.nx.centrality.eigenvector_centrality_numpy(a.G)
    except:
        eigCent = dict(zip(a.G.nodes(), ['NA' for n in a.G.nodes()]))
    extras.writeResults(eigCent, "eigCentNP", ofb, propDict=propDict, append=appVal)
    del(eigCent)
    
    eln = extras.edgeLengths(a.G, nodeWise=True)
    extras.writeResults(eln, "eln", ofb, propDict=propDict, append=appVal)
    del(eln)
    
    el = extras.edgeLengths(a.G)
    meanEL = np.mean(np.array((el.values()), dtype=float))
    extras.writeResults(meanEL, "meanEL", ofb, propDict=propDict, append=appVal)

    medianEL = np.median(np.array((el.values()), dtype=float))
    extras.writeResults(medianEL, "medianEL", ofb, propDict=propDict, append=appVal)

    del(el, meanEL, medianEL)
    
    # modularity metrics
    ci = bct.modularity_louvain_und(a.bctmat)
    Q = ci[1]
    ciN = a.assignbctResult(ci[0])
    extras.writeResults(Q, "Q", ofb, propDict=propDict, append=appVal)
    extras.writeResults(ciN, "ci", ofb , propDict=propDict, append=appVal)  
    
    pcCent = bct.participation_coef(a.bctmat,ci[0])
    pcCent = a.assignbctResult(pcCent)
    extras.writeResults(pcCent, "pcCent", ofb, propDict=propDict, append=appVal)
    del pcCent
    
    wmd = extras.withinModuleDegree(a.G, ciN)
    extras.writeResults(wmd, "wmd", ofb, append=appVal)
    del wmd    
    
    nM = len(np.unique(ci[0]))
    extras.writeResults(nM, "nM", ofb, propDict=propDict, append=appVal)
    del(nM)
    del(ci, ciN, Q)
    
    # rich club measures
    rc = mbt.nx.rich_club_coefficient(a.G, normalized=False)
    extras.writeResults(rc, "rcCoeff", ofb, propDict=propDict, append=appVal)
    del(rc)
    # robustness
    rbt = a.robustness()
    extras.writeResults(rbt, "robustness", ofb, propDict=propDict, append=False)

    # append any further iterations
    appVal = True
#    # infomap partitioning
#    bIM = infomap.nx2infomap(a.G)
#    del(bIM)
#    try:
#        f = open("nxdigraph.clu", "r") # recapture results from output file
#        modules = mbt.np.array([int(v.strip('\n')) for v in f.readlines()[1:]])
#        f.close()
#        remove("nxdigraph.clu")
#
#        ciNIM = a.assignbctResult(modules)
#        QIM = community.modularity(ciNIM, a.G)            
#        
#        pcCentIM = bct.participation_coef_sign(a.bctmat, modules)
#        pcCentIM = a.assignbctResult(pcCentIM)
#        wmdIM = extras.withinModuleDegree(a.G, ciNIM, weight='weight')
#        nMIM = len(np.unique(ciNIM[0]))
#    
#    except:
#        modules = np.array((np.repeat(np.nan, len(a.G.nodes()))))
#        QIM = "NA"
#        pcCentIM = {v:"NA" for v in a.G.nodes()}
#        wmdIM = {v:"NA" for v in a.G.nodes()}
#        ciNIM = {v:"NA" for v in a.G.nodes()}
#        nMIM = 0
#    
#    extras.writeResults(QIM, "QIM", ofbT, propDict=propDict, append=appValT)
#    extras.writeResults(ciNIM, "ciIM", ofbT, propDict=propDict, append=appValT)
#    del QIM
#    
#    extras.writeResults(nMIM, "nMIM", ofbT, propDict=propDict, append=appValT)
#    del nMIM
#    
#    extras.writeResults(pcCentIM, "pcCentIM", ofbT, propDict=propDict, append=appValT)
#    del pcCentIM
#
#    extras.writeResults(wmdIM, "wmdIM", ofbT, propDict=propDict, append=appValT)
#    del(wmdIM, ciNIM)


#    # Newman partitioning
#    try:
#        ciNm = bct.modularity_und(a.bctmat)
#        QNm= ciNm[1]
#        ciNNm = a.assignbctResult(ciNm[0])
#        nMNm = len(np.unique(ciNm[0]))
#        pcCentNm = bct.participation_coef_sign(a.bctmat,ciNm[0])
#        pcCentNm = a.assignbctResult(pcCentNm)
#        wmdNm = extras.withinModuleDegree(a.G, ciNNm, weight='weight')
#    
#    except:
#        ciNm=None
#        QNm = "NA"
#        ciNNm = dict(zip(a.G.nodes(), ['NA' for n in a.G.nodes()]))
#        nMNm = "NA" 
#        pcCentNm = dict(zip(a.G.nodes(), ['NA' for n in a.G.nodes()]))
#        wmdNm = dict(zip(a.G.nodes(), ['NA' for n in a.G.nodes()]))
#    
#    extras.writeResults(QNm, "QNm", ofbT, propDict=propDict, append=appValT)
#    extras.writeResults(ciNNm, "ciNm", ofbT, propDict=propDict, append=appValT)
#    del QNm
#    
#    extras.writeResults(nMNm, "nMNm", ofbT, propDict=propDict, append=appValT)
#    del nMNm
#    extras.writeResults(pcCentNm, "pcCentNm", ofbT, propDict=propDict, append=appValT)
#    
#    del pcCentNm
#    
#    extras.writeResults(wmdNm, "wmdNm", ofbT, propDict=propDict, append=appValT)
#    del(wmdNm,ciNNm,ciNm)

# weighted measures
a.importAdjFile(adjMatFile, delimiter=delim, excludedNodes=excludedNodes)
a.applyThreshold()
a.removeUnconnectedNodes()

##a.adjMatThresholding(MST=False)
a.weightToDistance()
ofb = '_'.join(["brain", "d"+dVal+"_"])
appVal = False
a.importSpatialInfo(parcelFile)  # read spatial information
a.makebctmat()

# weighted hub metrics
degs = a.G.degree(weight='weight')
extras.writeResults(degs, "degree_wt", ofb, append=appVal)

betCent = mbt.nx.centrality.betweenness_centrality(a.G, weight='distance')
extras.writeResults(betCent, "betCent_wt", ofb, append=appVal)

closeCent = mbt.nx.centrality.closeness_centrality(a.G, distance='distance')
extras.writeResults(closeCent, "closeCent_wt", ofb, append=appVal)

#hs = extras.hubscore(a.G, bc=betCent, cc=closeCent, degs=degs, weighted=True)
#extras.writeResults(hs, "hubscores_wt", ofb, append=appVal)
#del(hs, betCent)

eigCent = mbt.nx.centrality.eigenvector_centrality_numpy(a.G)
extras.writeResults(eigCent, "eigCentNP_wt", ofb, append=appVal)
del(eigCent)

# weighted modularity metrics
ci = bct.modularity_louvain_und_sign(a.bctmat)
Q = ci[1]
ciN = a.assignbctResult(ci[0])
extras.writeResults(Q, "Q_wt", ofb, append=appVal)
extras.writeResults(ciN, "ci_wt", ofb, append=appVal)  

nM = len(np.unique(ci[0]))
extras.writeResults(nM, "nM_wt", ofb, append=appVal)
del(nM)

wmd = extras.withinModuleDegree(a.G, ciN, weight='weight')
extras.writeResults(wmd, "wmd", ofb, append=appVal)
del wmd

clustCoeff = mbt.nx.average_clustering(a.G, weight="weight")
extras.writeResults(clustCoeff, "clusterCoeff_wt", ofb, append=appVal)
del(clustCoeff)

pl = mbt.nx.average_shortest_path_length(a.G, weight="distance")
extras.writeResults(pl, "pl_wt", ofb, append=appVal)
del(pl)

ge = extras.globalefficiency(a.G, weight="distance")
extras.writeResults(ge, "ge_wt", ofb, append=appVal)
del(ge)

le = extras.localefficiency(a.G, weight="distance")
extras.writeResults(le, "leN_wt", ofb, append=appVal)
del(le)

pcCent = np.zeros((len(a.G.nodes()), 10))
betCentT = np.zeros((len(a.G.nodes()), 10))
nM = np.zeros((10))
wmd = np.zeros((len(a.G.nodes()), 10))
Q = np.zeros((10))

appValT=False

#pcCentIM = np.zeros((len(a.G.nodes()), 10))
#nMIM = np.zeros((10))
#wmdIM = np.zeros((len(a.G.nodes()), 10))
#QIM = np.zeros((10))
#
#pcCentNm = np.zeros((len(a.G.nodes()), 10))
#nMNm = np.zeros((10))
#wmdNm = np.zeros((len(a.G.nodes()), 10))    
#QNm = np.zeros((10))


for n,i in enumerate([v for v in range(1,11)]):
    a.localThresholding(edgePC=i)
    a.removeUnconnectedNodes()
    a.makebctmat()
    a.weightToDistance()
    ofbT = '_'.join(["brain", thresholdtype, str(i), "d"+dVal+"_"])
    propDict = {"edgePC":a.edgePC}
    
    # weighted modularity metrics
    ci = bct.modularity_louvain_und_sign(a.bctmat)
    QWt = ci[1]
    Q[n] = QWt
    ciN = a.assignbctResult(ci[0])
    extras.writeResults(QWt, "QWt", ofbT, propDict=propDict, append=appValT)
    extras.writeResults(ciN, "ciWt", ofbT, propDict=propDict, append=appValT)
    del QWt
    
    nMWt = len(np.unique(ci[0]))
    nM[n] = nMWt
    extras.writeResults(nMWt, "nMWt", ofbT, propDict=propDict, append=appValT)
    del(nMWt)

    wmdWt = extras.withinModuleDegree(a.G, ciN, weight='weight')
    wmd[:,n] = [wmdWt[v] for v in a.G.nodes()]
    extras.writeResults(wmdWt, "wmdWt", ofbT, propDict=propDict, append=appValT)
    del wmdWt
    
    pcCentWt = bct.participation_coef_sign(a.bctmat,ci[0])
    pcCent[:,n] = pcCentWt
    pcCentWt = a.assignbctResult(pcCentWt)
    extras.writeResults(pcCentWt, "pcCentWt", ofbT, propDict=propDict, append=appValT)
    
    bcT = mbt.nx.centrality.betweenness_centrality(a.G, weight='distance')    
    betCentT[:,n] = [bcT[v] for v in a.G.nodes()]
    appValT=True
    
#    # infomap partitioning
#    bIM = infomap.nx2infomap(a.G)
#    del(bIM)
#    f = open("nxdigraph.clu", "r") # recapture results from output file
#    modules = mbt.np.array([int(v.strip('\n')) for v in f.readlines()[1:]])
#    f.close()
#    remove("nxdigraph.clu")
#    
#    ciNIM = a.assignbctResult(modules)
#    QIMWt = community.modularity(ciNIM, a.G)
#    QIM[n] = QIMWt
#    extras.writeResults(QIMWt, "QIMWt", ofbT, propDict=propDict,append=appValT)
#    extras.writeResults(ciNIM, "ciIMWt", ofbT, propDict=propDict, append=appValT)
#    del(QIMWt)
#    
#    nMIMWt = len(np.unique(modules))
#    nMIM[n] = nMIMWt
#    extras.writeResults(nMIMWt, "nMIMWt", ofbT, propDict=propDict, append=appValT)
#    del(nMIMWt)
#    
#    pcCentIMWt = bct.participation_coef_sign(a.bctmat, modules)
#    pcCentIM[:,n] = pcCentIMWt
#    pcCentIMWt = a.assignbctResult(pcCentIMWt)
#    extras.writeResults(pcCentIMWt, "pcCentIMWt", ofbT, propDict=propDict, append=appValT)
#    del(pcCentIMWt)
#    
#    wmdIMWt = extras.withinModuleDegree(a.G, ciNIM, weight='weight')
#    wmdIM[:,n] = [wmdIMWt[v] for v in a.G.nodes()]
#    extras.writeResults(wmdIMWt, "wmdIMWt", ofbT, propDict=propDict, append=appValT)
#    del wmdIMWt

#    # Newman partitioning
#    ciNm = bct.modularity_und(a.bctmat)
#    QNmWt = ciNm[1]
#    QNm[n] = QNmWt
#    ciNNm = a.assignbctResult(ciNm[0])
#    extras.writeResults(QNmWt, "QNmWt", ofbT, propDict=propDict, append=appValT)
#    extras.writeResults(ciNNm, "ciNmWt", ofbT, propDict=propDict, append=appValT)  
#    
#    nMNmWt = len(np.unique(ciNm[0]))
#    nMNm[n] = nMNmWt
#    extras.writeResults(nMNmWt, "nMNmWt", ofbT, propDict=propDict, append=appValT)
#    del(nMNmWt)
#
#    pcCentNmWt = bct.participation_coef_sign(a.bctmat,ciNm[0])
#    pcCentNm[:,n] = pcCentNmWt
#    pcCentNmWt = a.assignbctResult(pcCentNmWt)
#    extras.writeResults(pcCentNmWt, "pcCentNmWt", ofbT, propDict=propDict, append=appValT)
#    
#    wmdNmWt = extras.withinModuleDegree(a.G, ciNNm, weight='weight')
#    wmdNm[:,n] = [wmdNmWt[v] for v in a.G.nodes()]
#    extras.writeResults(wmdNmWt, "wmdNmWt", ofbT, propDict=propDict, append=appValT)
#    del wmdNmWt


Q = np.mean(Q)
extras.writeResults(Q, "QWt_wt", ofb, append=appVal)
del(Q)

pcCent = a.assignbctResult(np.mean(pcCent, axis=1))
extras.writeResults(pcCent, "pcCentWt_wt", ofb, append=appVal)
del(pcCent,ci)

betCentT = a.assignbctResult(np.mean(betCentT, axis=1))
extras.writeResults(betCentT, "betCentWtT_wt", ofb, append=appVal)

nM = np.mean(nM)
extras.writeResults(nM, "nMWt_wt", ofb, append=appVal)
del(nM)

wmd = a.assignbctResult(np.mean(wmd, axis=1))
extras.writeResults(wmd, "wmdWt_wt", ofb, append=appVal)
del(wmd)


#hs = extras.hubscore(a.G, bc=betCentT, cc=closeCent, degs=degs, weighted=True)
#extras.writeResults(hs, "hsT_wt", ofb, append=appVal)
#del(hs, betCentT)

## Infomap
#QIM = np.mean(QIM)
#extras.writeResults(QIM, "QIMWt_wt", ofb, append=appVal)
#del(QIM)
#
#pcCentIM = a.assignbctResult(np.mean(pcCentIM, axis=1))
#extras.writeResults(pcCentIM, "pcCentIMWt_wt", ofb, append=appVal)
#del(pcCentIM)
# 
#wmdIM = a.assignbctResult(np.mean(wmdIM, axis=1))
#extras.writeResults(wmdIM, "wmdIMWt_wt", ofb, append=appVal)
#del(wmdIM)
#
#nMIM = np.mean(nMIM)
#extras.writeResults(nMIM, "nMIMWt_wt", ofb, append=appVal)
#del(nMIM)

## Newman
#QNm = np.mean(QNm)
#extras.writeResults(QNm, "QNmWt_wt", ofb, append=appVal)
#del(QNm)
#
#pcCentNm = a.assignbctResult(np.mean(pcCentNm, axis=1))
#extras.writeResults(pcCentNm, "pcCentNmWt_wt", ofb, append=appVal)
#del(pcCentNm)
# 
#wmdNm = a.assignbctResult(np.mean(wmdNm, axis=1))
#extras.writeResults(wmdNm, "wmdNmWt_wt", ofb, append=appVal)
#del(wmdNm)
#
#nMNm = np.mean(nMNm)
#extras.writeResults(nMNm, "nMNmWt_wt", ofb, append=appVal)
#del(nMNm)

