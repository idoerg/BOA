"""
Code taken from pydoc.net

Big note: The ROC curves expect a probability vector for the predicted values
"""
import rpy2.robjects as robjects

robjects.r("library(ROCR)")

def device(filename):
    robjects.r.pdf(file=filename,width=5,height=5)

def close():
    robjects.r("dev.off()")

def transform(labels,refLabel):
    return [x==refLabel for x in labels]
def transformDistribution(pds,refLabel):
    return [ x.prob(refLabel) for x in pds if x!=None]
def plot(ref_labels,pred_labels,title,y="tpr",x="fpr",diag=True):
    try:
        reference = robjects.IntVector(ref_labels)
    except:
        print "Problem with reference data"
    else:
        try:
            predictions= robjects.FloatVector(pred_labels)
        except:
            print "Problem with predictive data"
        else:
            pred = robjects.r.prediction(predictions,reference)
            perf = robjects.r.performance(pred,y,x)
            auc  = robjects.r.performance(pred,"auc").do_slot("y.values")[0][0]
            print "AUC",auc
            #title="%s (AUC=%f)"%(title,auc)
            try:
                robjects.r.plot(perf,main=title,colorize=True)
            except:
                print "problem with plotting"
            if diag:
                robjects.r.abline(a=0,b=1)
