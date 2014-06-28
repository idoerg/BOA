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
import rforests
import textseq_sim

root = os.environ['BACFINDER_HOME']
trainDir = "%s/data/training/protein"%root
labelFile = "%s/data/training/training_proteins4.txt"%root
zip = "test_serial.zip"
function_check="function_check.txt"
db="sim.db"
fasta="sim.fa"
"""
nb = nbayes.NBayes(trainDir,labelFile)
nb.train()
ref_labels,pred_labels = nb.testClassify(70)
ref_labels  = [x for x in ref_labels if x!='na']
pred_labels  = [x for x in pred_labels if x!='na']


reft_labels = roc.transform(ref_labels,'null')
null_labels = roc.transformDistribution(pred_labels,'null')
roc.device("roc.pdf")
roc.plot(reft_labels,null_labels,'Null Class Naive Bayes ROC curve')
reft_labels = roc.transform(ref_labels,'toxin')
toxin_labels = roc.transformDistribution(pred_labels,'toxin')
roc.plot(reft_labels,toxin_labels,'Toxin Class Naive Bayes ROC curve')
reft_labels = roc.transform(ref_labels,'modifier')
modifier_labels = roc.transformDistribution(pred_labels,'modifier')
roc.plot(reft_labels,modifier_labels,'Modifier Class Naive Bayes ROC curve')
reft_labels = roc.transform(ref_labels,'transport')
transport_labels = roc.transformDistribution(pred_labels,'transport')
roc.plot(reft_labels,transport_labels,'Transport Class Naive Bayes ROC curve')
reft_labels = roc.transform(ref_labels,'immunity')
immunity_labels = roc.transformDistribution(pred_labels,'immunity')
roc.plot(reft_labels,immunity_labels,'Immunity Class Naive Bayes ROC curve')
reft_labels = roc.transform(ref_labels,'regulator')
regulator_labels = roc.transformDistribution(pred_labels,'regulator')
roc.plot(reft_labels,regulator_labels,'Regulator Class Naive Bayes ROC curve')

mnb = rforests.RForests(trainDir,labelFile,numTrees=500)
ref_labels,pred_labels = mnb.testClassify(70)
ref_labels  = [x for x in ref_labels if x!='na']
reft_labels = roc.transform(ref_labels,'null')
null_labels = roc.transformDistribution(pred_labels,'null')
roc.plot(reft_labels,null_labels,'Null Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'toxin')
toxin_labels = roc.transformDistribution(pred_labels,'toxin')
roc.plot(reft_labels,toxin_labels,'Toxin Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'modifier')
modifier_labels = roc.transformDistribution(pred_labels,'modifier')
roc.plot(reft_labels,modifier_labels,'Modifier Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'transport')
transport_labels = roc.transformDistribution(pred_labels,'transport')
roc.plot(reft_labels,transport_labels,'Transport Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'immunity')
immunity_labels = roc.transformDistribution(pred_labels,'immunity')
roc.plot(reft_labels,immunity_labels,'Immunity Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'regulator')
regulator_labels = roc.transformDistribution(pred_labels,'regulator')
roc.plot(reft_labels,regulator_labels,'Regulator Class Random Forests ROC curve')
"""
#accs = nb.learningCurve(numTrials=10)

#plt.plot(xrange(len(accs)),accs)
#plt.xlabel("Training set size")
#plt.ylabel("Cross-Validation Accuracy")
#plt.title("Learning Curve")
#plt.show()

#Create ROC curve with simulation data
"""
if os.path.exists(db):
    os.remove(db)
if os.path.exists(fasta):
    os.remove(fasta)
    
nsim = textseq_sim.TextSeqSim(trainDir,labelFile)
nsim.simulationOutput(10000,function_check,db,fasta)

mnb = rforests.RForests(trainDir,labelFile,numTrees=500)
mnb.train()
ref = nsim.readOutput(function_check)

prot_ids,ref_labels = ref
prot_ids,pred_labels,pds = mnb.prob_classify(db,fasta)
ref_labels  = [x for x in ref_labels if x!='na' and x!='regulator' and x!='transport']

reft_labels = roc.transform(ref_labels,'null')
null_labels = roc.transformDistribution(pds,'null')
roc.device("roc.pdf")
roc.plot(reft_labels,null_labels,'Null Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'toxin')
toxin_labels = roc.transformDistribution(pds,'toxin')
roc.plot(reft_labels,toxin_labels,'Toxin Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'modifier')
modifier_labels = roc.transformDistribution(pds,'modifier')
roc.plot(reft_labels,modifier_labels,'Modifier Class Random Forests ROC curve')
reft_labels = roc.transform(ref_labels,'immunity')
immunity_labels = roc.transformDistribution(pds,'immunity')
roc.plot(reft_labels,immunity_labels,'Immunity Class Random Forests ROC curve')
"""
rf = rforests.RForests(trainDir,labelFile,numTrees=500)
rf.train()
accs = rf.learningCurve(numTrials=10)
print accs
plt.plot(xrange(len(accs)),accs)
plt.xlabel("Training set size")
plt.ylabel("Cross-Validation Accuracy")
plt.title("Learning Curve")
plt.show()

#roc.close()



