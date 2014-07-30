"""
Perform expectation maximization algorithm to figure out the mixture model that best represents the data
The data is made up of true bacteriocin associated genes and not bacteriocin associated genes.
Each data point is a score provided by HMMER.  Hence, the true and false distributions can be
represented as a mixture model

It is hypothesized that the mixture model is composed of a Gumbel distribution and another distribution
For the sake of simplicity, lets first model this as a Gaussian Mixture Model.
"""
import numpy
from scipy.stats import gumbel_r
from scipy.stats import norm
from numpy import random
import os,site,sys
from math import *

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import hmmer

class GaussianMixtureModel(object):
	"""
	mu1 and mu2 are mean variables
	sig2 is the standard deviation for the normal distribution
	sig1 is the scaling parameter for the gumbel distribution
	"""
	def __init__(self,scores=None):
		self.scores = scores
	""" Read in HMMER scores """
	def parse(hmmin):
		self.scores = hmmer.parse_scores(hmmin)

	def pdf_model(self,x,params,funcs=[norm.pdf,norm.pdf]): 
		assert len(funcs)==2 #Only can have two functions
		assert len(params)==5 
		mu1, sig1, mu2, sig2, pi_1 = params 
		return pi_1*funcs[0](x,mu1,sig1)+(1-pi_1)*funcs[1](x,mu2,sig2)
	""" This code was stolen from the following website 
	http://nbviewer.ipython.org/github/tritemio/notebooks/blob/master/Mixture_Model_Fitting.ipynb
	"""
	def expectation_maximization(self,max_iter=1000):
		# Initial guess of parameters and initializations
		p0 = numpy.array([-0.2,0.2,0.8,0.2,0.5])

		mu1, sig1, mu2, sig2, pi_1 = p0
		mu = numpy.array([mu1, mu2])
		sig = numpy.array([sig1, sig2])
		pi_ = numpy.array([pi_1, 1-pi_1])
		s = self.scores
		gamma = numpy.zeros((2, s.size))
		N_ = numpy.zeros(2)
		p_new = p0
		func = [norm.pdf,norm.pdf] 
		# EM loop
		counter = 0
		converged = False
		N = len(self.scores)
		while not converged:
		    # Compute the responsibility func. and new parameters
		    for k in [0,1]:
		    	pdf_model = self.pdf_model(s, p_new,func) 
		        gamma[k,:] = pi_[k]*func[k](s, mu[k], sig[k])/ pdf_model
		        N_[k] = 1.*gamma[k].sum()
		        mu[k] = sum(gamma[k]*s)/N_[k]
		        sig[k] = sqrt( sum(gamma[k]*(s-mu[k])**2)/N_[k] )
		        pi_[k] = N_[k]/s.size
		    p_new = [mu[0], sig[0], mu[1], sig[1], pi_[0]]
		    assert abs(N_.sum() - N)/float(N) < 1e-6 
		    assert abs(pi_.sum() - 1) < 1e-6
		    
		    # Convergence check
		    counter += 1
		    converged = counter >= max_iter
		print "Means:   %6.3f  %6.3f" % (p_new[0], p_new[2])
		print "Std dev: %6.3f  %6.3f" % (p_new[1], p_new[3])
		print "Mix (1): %6.3f " % p_new[4]
		return p_new

if __name__=="__main__":
	import unittest
	class TestMixtureModel(unittest.TestCase):
		def setUp(self):
			N = 1000
			a = 0.3
			s1 =  random.normal(0, 0.08, size=N*a)
			s2 = random.normal(0.6,0.12, size=N*(1-a))
			s = numpy.concatenate([s1,s2])
			self.data = s
		def test1(self):
			gmm = GaussianMixtureModel(self.data)
			params = gmm.expectation_maximization(1000)
			self.assertTrue(abs(params[0])<0.01)
			self.assertTrue(abs(params[1]-0.08)<0.01)
			self.assertTrue(abs(params[2]-0.6)<0.01)
			self.assertTrue(abs(params[3]-0.12)<0.01)
			self.assertTrue(abs(params[4]-0.3)<0.01)

	unittest.main()


