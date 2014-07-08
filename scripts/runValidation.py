




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

if __name__=="__main__":
    
    operons = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/intermediate/operons.txt"
    ptt = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/data/all.ptt"
    
    categorize(operons,ptt,"categories.txt") 