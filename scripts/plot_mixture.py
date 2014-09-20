"""
Plots gaussian mixture model
"""
import re

import numpy
from scipy.stats import gumbel_r
from scipy.stats import norm
from numpy import random
import os,site,sys
from math import *
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import hmmer
import pylab
from mixture import GaussianMixtureModel
import matplotlib.pyplot as plt
from sklearn import mixture


def filter_bagel(bagel,Gbagel):
    scores = []
    genome_bagels = set()
    with open(Gbagel,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = re.split("\s+",ln)
            _id = toks[0]
            genome_bagels.add(_id)
    with open(bagel,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = re.split("\s+",ln)
            _id = toks[2]
            if _id in genome_bagels:
                scores.append(float(toks[5]))
    return scores



#folder = "/Users/mortonyt/Documents/MiamiBio/workspace"
folder = "/media/HD/Documents/Jamie/MiamiBio/Bacfinder/workspace/gmm_workspace"
toxin_scores     = hmmer.parse_scores("%s/toxin_mod.out"%folder)
#modifier_scores  = hmmer.parse_scores("%s/modifier.out"%folder)
#immunity_scores  = hmmer.parse_scores("%s/immunity.out"%folder)
#regulator_scores = hmmer.parse_scores("%s/regulator.out"%folder)
#transport_scores = hmmer.parse_scores("%s/transport.out"%folder)
bagel_scores     = filter_bagel("%s/bagel_reduced.out"%folder,
                                "%s/bagel_whole_genome.out"%folder)
#all_scores = numpy.array(toxin_scores+modifier_scores+immunity_scores+regulator_scores+transport_scores)
all_scores = toxin_scores
print all_scores[:10]
print bagel_scores[:10]
print "Number of all scores",len(all_scores)
print "Number of bagel scores",len(bagel_scores)
print "Min bagel score",min(bagel_scores)
# gmm = GaussianMixtureModel(all_scores)
# g = mixture.GMM(n_components=2)
# model = g.fit(all_scores)
# x = numpy.linspace(min(all_scores),max(all_scores),1000)
# logprob,reps = model.score_samples(x)
# pdf = numpy.exp(logprob)
# indiv_pdf = reps*pdf[:,numpy.newaxis]

print "Calculated pdf"
n, bins, patches=plt.hist(all_scores, 100, normed=True,histtype="stepfilled")
hmmer=plt.setp(patches, 'facecolor', 'r', 'alpha', 0.5,)
n, bins, patches=plt.hist(bagel_scores, 30,normed=True,histtype="stepfilled")
bagel=plt.setp(patches, 'facecolor', 'b', 'alpha', 0.5)
plt.legend(['HMMER scores','Bagel scores'])
plt.xlabel('Score')
plt.ylabel('Normalized Counts')
plt.xlim(0,1000)
plt.ylim(0,0.02)
#plt.plot(x, pdf, '-k')
#plt.plot(x, indiv_pdf, '--r')
plt.show()



