"""
Plots ROC curves for classifiers 
"""
import os,site,sys

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cbook as cbook
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage 
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import roc
import dtree 
import nbayes
import mnbayes

root = os.environ['BACFINDER_HOME']
trainDir = "%s/data/training/protein"%root
labelFile = "%s/data/training/training_proteins.txt"%root
zip = "test_serial.zip"
nb = nbayes.NBayes(trainDir,labelFile)
nb.train()

ref_labels,pred_labels = nb.testClassify(50)
ref_labels  = [x for x in ref_labels if x!='na']

ref_labels = roc.transform(ref_labels,'null')
pred_labels = roc.transformDistribution(pred_labels,'null')


roc.device("null_nbayes_roc.pdf")
roc.plot(ref_labels,pred_labels,'Null Class ROC curve')

"""
dt = dtree.DTree(trainDir,labelFile)
dt.train()

ref_labels,pred_labels = dt.testClassify(50)
ref_labels  = [x for x in ref_labels if x!='na']
pred_labels = [x for x in pred_labels if x!='na']

ref_labels = roc.transform(ref_labels,'null')
pred_labels = roc.transform(pred_labels,'null')
roc.device("null_dtree_roc.pdf")
roc.plot(ref_labels,pred_labels,'Null Class ROC curve')
"""

#mnb = mnbayes.MNBayes(trainDir,labelFile)
#mnb.rocCurve()

#mnb.train()
#ref_labels,pred_labels = mnb.testClassify(50)
#ref_labels  = [x for x in ref_labels if x!='na']
#ref_labels = roc.transform(ref_labels,'null')
#pred_labels = roc.transformDistribution(pred_labels,'null')
#roc.device("null_mnbayes_roc.pdf")
#roc.plot(ref_labels,pred_labels,'Null Class ROC curve')


#accs = nb.learningCurve(numTrials=10)

#plt.plot(xrange(len(accs)),accs)
#plt.xlabel("Training set size")
#plt.ylabel("Cross-Validation Accuracy")
#plt.title("Learning Curve")
#plt.show()
roc.close()



