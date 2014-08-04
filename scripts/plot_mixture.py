"""
Plots gaussian mixture model
"""


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


def best_scores(hits):
	prev = hits[0]
	prev_name,_ = prev
	buf = []
	maxScores = []
	for hit in hits:
		name,score = hit
		if name==prev_name:
			buf.append(score)
		else:
			maxScores.append(max(buf))
			buf = [score]
			prev_name = name
	return maxScores

folder = "/Users/mortonyt/Documents/MiamiBio/workspace"
boa_scores     = hmmer.parse_scores("%s/boa_scores.out"%folder)
bagel_scores     = hmmer.parse_scores("%s/bagel_scores.out"%folder)
boa_scores = best_scores(boa_scores)
bagel_scores = best_scores(bagel_scores)

print "Number of boa scores",len(boa_scores)
print "Number of bagel scores",len(bagel_scores)

gmm = GaussianMixtureModel(boa_scores)
#params = gmm.expectation_maximization(1000)
g = mixture.GMM(n_components=2)
model = g.fit(boa_scores)
x = numpy.linspace(min(boa_scores),max(boa_scores),1000)
logprob,reps = model.score_samples(x)
pdf = numpy.exp(logprob)
indiv_pdf = reps*pdf[:,numpy.newaxis]
print "Calculated pdf"
n, bins, patches=plt.hist(boa_scores, 30, normed=True,histtype="stepfilled")
plt.setp(patches, 'facecolor', 'c', 'alpha', 0.5)
n, bins, patches=plt.hist(bagel_scores, 30, normed=True,histtype="stepfilled")
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.5)

#plt.plot(x, pdf, '-k')
#plt.plot(x, indiv_pdf, '--r')
print "Display"
plt.show()



