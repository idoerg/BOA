"""
This module conducts a majority vote on each cluster to determine the function
There should only be 7 possible function classes: toxins, modifiers, regulators, transporters, immunity, null, NA (not assigned)
"""
import os,site,sys
import re
from collections import Counter
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import cdhit
cluster_reg = re.compile(r">(\S+)")

class MVote(object):
    def __init__(self,titles,labels,clusters):
        self.titles = titles
        self.labels = labels
        self.clusters = clusters
        self.labeldict = {self.titles[i]:self.labels[i] 
                          for i in range(len(self.titles))}
    """ Assign a function to a cluster number """
    def assign(self):
        clusterNum = 0
        labels = []
        for cluster in self.clusters:
            function = self.vote(cluster)
            labels.append( (clusterNum,function) )
            clusterNum+=1
        return labels
    """ Determine the most commonly seen function"""
    def vote(self,cluster):
        funcCounts = Counter()
        assigned = 0
        for seq_id in cluster:
            title = cluster_reg.findall(seq_id)[0]
            toks = title.split('|')
            protid = toks[-1]
            protid = protid[:-3] #Remove shitty cdhit labels
            print "Seq id",seq_id
            print "Protein id",protid
            function = self.labeldict[protid]
            if function!='na': #Don't count unassigned functions
                funcCounts[function]+=1
                assigned+=1
        if assigned==0: return 'na'
        func,cnts = funcCounts.most_common(1)[0]
        return func
        
    
if __name__=="__main__":
    import unittest
    
    class TestMajorityVote(unittest.TestCase):
        def setUp(self):
            self.testfa = "test.fa"
            self.testclr = "testclr.fa"
            self.seqs = [
                    '>ALQPTKWCSI|_|_|_|_|YRQNRAKNBT.1',
                    'FAVMRYCQSVADTEWFNQQWENMSYAWQSTRRVLRDEWEQMEHEWDTLSTWMECANVNHHVGLMLIKQNGVTLYPVPFMPCKWIKEMVMCWTEQFQTNRLMQMQTTLLPQYAKEQWSWPQEHKNLAIEEWPECVDTYTFWMDGMINERQRTVAAAIAWGLNFRPQYKSVDTPVFHGTAFKCADDMKHFINWLRIGNFSNT',
                    '>GLMMCOTXII|_|_|_|_|XSDJNGQEAJ.1',
                    'KRIGKPWFTGVWHYLFKSQPIENMVIWENEQVWRQERTSNIMWPPTDTAYIQWLAQPCGPNLEQLYKNLKNVYYWIPNCANMWSMMGIWSQMGIAGGIMVRDSKWRVLYQPYCIGSLQPWNEDYLHPVSSKGLHLMIRYGHCYHQSHKAEPQDDCGFMDLYERWESGMGRCPCGTWDALDRTSNASMDPGNNQKILVPSF',      
                    '>HASTETYBNW|_|_|_|_|BHSJYDJDLO.1',
                    'DQQSKRSIRDKYMGRMWYRCTWHLYVHYKPDAAWILCALVAVYATTPSHLVQRYWAQTDRHDNGHGGKVDLQHTESPQVLIWHSIQADMMKIVFLVFMIMFVRVYIKCFPAYIPRCCLGECGDGTNIMAFFFYDAWTSNTIMQDLKQLFDQWVITMGSTSRSWYDCLAWSLWMHMQQWHFCGFGPRTTHHSEAIWLKDGE',
                    '>HBDVQKSYLI|_|_|_|_|KGLKXRIFIU.1',
                    'EVPMCMGTIRTTGRNMGGAVRGYSVRTMHEDVEQKSEAPPTTLKITSPMYYMFNPKPTMKWIWGIEMFDHFRRKRQIDIYPAYWGGTPQTQEANCKNTDPQAICLHEIWHLPNQLCGNQPPAIAINTVMCVRYCEMWGRDRGKPGEYPPCHTPHEMEHYLWERQSTYECETFFYVKFTGTNPCFGAWTHQVNGDYMIMKD',
                    '>XKIXHJDPDO|_|_|_|_|CUEJZIKBRC.1',
                    'MNLDTWNKWPVPGIQHIARVAMNLGRHNCWSLMFCIECWGVTCSEAWWESNYHCELRSQCTPCAKFWVPNGNCARAYICDLGHLSLHWRCQAWSYCHGNESRQNMSQFQCMVPPLQQDMQPRLFWRLPWMNESCWFHVYKINVYTTMLMVIMTDLSDHEDKWFLSPTDPIMQGPENEFDIHAVEQKSKFVCPSQKMNPAQ',
                    '>XHSLIXVLRY|_|_|_|_|UMNKQPFGQW.1',
                    'MRVQYDEANCNSYGCQVGVLKTRNYMSATELRAYVQLNVWDWSRWRRIFDIESQTCGDNAMRGYLNDKSARPFMWSWYLVYIKPLVEWLPYDTWISCDSNSMHPNSFSWRFGPNQSMMPMALAFLGGVDCSTRMNEQLIWWAAQTDDFPCTGHQYLSAIENNRNHSDIHPESCNCVHTNFEHWEKEYNCPKWPISKEYRP',
                    '>KRCBQLSAHX|_|_|_|_|YQQAUGVMLS.1',
                    'CWHSWLWMYIGLSRIVMVWYSISISSCRKLKLLDLYQWHVGARGECAYTPWEDGFVGEHFWFHMQGCRTNARDSCTKDIQFRGHCGFVGKNVTERGDETVLHTAKPEHTWLLFLINVHDHKVKIEKFLNTYMDFYLSQYMGQVFPSSKWDYCGEEFPVVGLCMKFRNLGMYTGNIAQSSYWVCVAPVHSYGKPVYMHWDD',
                    '>QDPZTHDCEM|_|_|_|_|ORRFAHZJOM.1',
                    'CWHSWLWMYIGLSRIVMVWYSISISSCRKLKLLDLYQWHVGARGECAYTPWEDGFVGEHFWFHMQGCRTNARDSCTKDIQFRGHCGFVGKNVTERGDETVLHTAKPEHTWLLFLINVHDHKVKIEKFLNTYMDFYLSQYMGQVFPSSKWDYCGEEFPVVGLCMKFRNLGMYTGNIAQSSYWVCVAPVHSYGKPVYMHWDD',
                    '>KRCBQLSAHX|_|_|_|_|YQQAUGVMLX.1',
                    'CWHSWLWMYIGLSRIVMVWYSISISSCRKLKLLDLYQWHVGARGECAYTPWEDGFVGEHFWFHMQGCRTNARDSCTKDIQFRGHCGFVGKNVTERGDETVLHTAKPEHTWLLFLINVHDHKVKIEKFLNTYMDFYLSQYMGQVFPSSKWDYCGEEFPVVGLCMKFRNLGMYTGNIAQSSYWVCVAPVHSYGKPVYMHWDD',
                    '>QDPZTHDCEM|_|_|_|_|ORRFAHZJOX.1',
                    'CWHSWLWMYIGLSRIVMVWYSISISSCRKLKLLDLYQWHVGARGECAYTPWEDGFVGEHFWFHMQGCRTNARDSCTKDIQFRGHCGFVGKNVTERGDETVLHTAKPEHTWLLFLINVHDHKVKIEKFLNTYMDFYLSQYMGQVFPSSKWDYCGEEFPVVGLCMKFRNLGMYTGNIAQSSYWVCVAPVHSYGKPVYMHWDD',
                    ]
            open(self.testfa,'w').write('\n'.join(self.seqs))
            self.titles = ['YRQNRAKNBT.1',
                          'XSDJNGQEAJ.1',
                          'BHSJYDJDLO.1',
                          'KGLKXRIFIU.1',
                          'CUEJZIKBRC.1',
                          'UMNKQPFGQW.1',
                          'YQQAUGVMLS.1',
                          'ORRFAHZJOM.1',
                          'YQQAUGVMLX.1',
                          'ORRFAHZJOX.1']
            self.labels= ['na',
                          'modifier',
                          'toxin',
                          'modifier',
                          'toxin',
                          'immunity',
                          'immunity',
                          'transport',
                          'immunity',
                          'na',
                          ]
        def testmvote(self):
            clr = cdhit.CDHit(self.testfa,self.testclr,0.7)
            clr.run()
            clr.parseClusters()
            mv = MVote(self.titles,
                       self.labels,
                       clr.clusters)
            results = mv.assign()
            self.assertEquals(len(results),7)
            correctLabels = [ 'na',
                              'modifier',
                              'toxin',
                              'modifier',
                              'toxin',
                              'immunity',
                              'immunity']
            correct = zip(range(7),correctLabels)
            self.assertEquals(set(correct),set(results))
    unittest.main()
    
    
    
    
    
    
    
    