




import os
import sys
import site
import re
import numpy as np
import numpy.random
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))

from hmmer_validation import *
from bagel_validation import *

if __name__=="__main__":
    
    operons = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/intermediate/operons.txt"
    predicted_operons = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/intermediate/predicted_operons.txt"
    bagel = "/home/mortonjt/Projects/Bacfinder/bacteriocins/bagel.csv"
    
    ptt = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/data/all.ptt"
    #bagel,both,detect = bagel_compare_species(operons,bagel)
    #open("bagel.out",'w').write('\n'.join(map(str,sorted(list(bagel))))+'\n')
    #open("both.out",'w').write('\n'.join(map(str,sorted(list(both))))+'\n')
    #open("detect.out",'w').write('\n'.join(map(str,sorted(list(detect))))+'\n')

    categorize(operons,ptt,"categories.txt")
    categorize(predicted_operons,ptt,"predicted_categories.txt")
