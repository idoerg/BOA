
import os,site,sys
from collections import *
# from bx.intervals import *
import re
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
# from annotated_genes import AnnotationTree

    

""" Compare discovered genes to Bagel """
def bagel_compare_species(operon_file,bagel_csv):
    bagel_species = set()
    detected_species = set()
    with open(bagel_csv,'r') as bagel_handle:
        for ln in bagel_handle:
            toks = ln.split(',')
            #... skip ...
            if toks[1]=="ID":continue
            if toks[3]=="":continue
            text = toks[3]
            species = " ".join(text.split(' ')[:2])
            bagel_species.add(species)
    with open(operon_file,'r') as operon_handle:
         for ln in operon_handle:
             if ln[0]=="-":continue
             ln = ln.rstrip()
             toks = ln.split("|")
             text = toks[7]
             species = " ".join(text.split(' ')[:2])
             detected_species.add(species)
    bagel  =  bagel_species - detected_species
    detect =  detected_species - bagel_species
    both   =  bagel_species & detected_species
    return bagel,both,detect

if __name__=="__main__":
    import unittest
    class TestBagelComparison(unittest.TestCase):
        def setUp(self):
            bagel_entries = [
                ",,,,,,,,,,,,,,,,,,,,,,,,",
                ",ID,NAME,ORGANISM,Modsystem,LiteratureRef,NCBI,Uniprot,Sequence,Subclass,Leaderlength,DOI,,http://dx.doi.org/,,Uniprot,,NCBI,,http://www.ncbi.nlm.nih.gov/protein/,,http://www.uniprot.org/uniprot/,,Comment,",
                ",,,,,,,,,,,,,,,,,,,,,,,,",
                ",,,,,,,,,,,,,,,,,,,,,,,,",
                ",1.1,Anacyclamide_(AcyE),Helicobacter pylori 26695,CyaG,REF,,D2K7B5_9NOST,MTKKNIRPQQVAPVERETISTAKDQSGQVQAQSSVIWGSPVPFAGDDAE,Cyanobactin,,,10.1128/AEM.01061-09,,http://dx.doi.org/10.1128/AEM.01061-09,,D2K7B5_9NOST,,,,http://www.ncbi.nlm.nih.gov/protein/,,http://www.uniprot.org/uniprot/D2K7B5_9NOST,,",
                ",2.1,Ancovenin,Helicobacter pylori 26695,LanM,REF,,LANC_STRS6,CVQSCSFGPLTWSCDGNTK,lanthipeptide B,,,10.1016/S0040-4039(00)89174-1,,http://dx.doi.org/10.1016/S0040-4039(00)89174-1,,LANC_STRS6,,,,,,http://www.uniprot.org/uniprot/LANC_STRS6,,P38655.1",
                ",3.1,Astexin,Asticcacaulis excentricus,LasB,REF,,E8RMD3,MHTPIISETVQPKTAGLIVLGKASAETRGLSQGVEPDIGQTYFEESRINQD,Lasso peptide,,,10.1073/pnas.1208978109,,http://dx.doi.org/10.1073/pnas.1208978109,,E8RMD3,,,,,,http://www.uniprot.org/uniprot/E8RMD3,,",
                ",4.1,Avermipeptin,Streptomyces avermitilis,LanKC,REF,,Q825F5,MALLDLQTMESDEHTGGGGASTVSLLSCVSAASVLLCL,lanthipeptide C,20,,10.1002/cbic.201200118,,http://dx.doi.org/10.1002/cbic.201200118,,Q825F5,,,,,,http://www.uniprot.org/uniprot/Q825F5,,",
                ",5.1,BhtA1,Streptococcus ratti,LanM,REF,73486983,,MKEIQKAGLQEELSILMDDANNLEQLTAGIGTTVVNSTFSIVLGNKGYICTVTVECMRNCQ,lanthipeptide B,29,,10.1016/j.femsle.2005.09.003,,http://dx.doi.org/10.1016/j.femsle.2005.09.003,,,,73486983,,http://www.ncbi.nlm.nih.gov/protein/AAZ76603.1,,http://www.uniprot.org/uniprot/,,AAZ76603.1"]

            operons = [ "AE000511.1_3|transport.fa.cluster10.fa|1.6e-20|5|213|803726|804395|Helicobacter pylori 26695, complete genome  ",
                        "AE000511.1_3|transport.fa.cluster12.fa|1.9e-34|15|217|803732|804374|Helicobacter pylori 26695, complete genome ",
                        "AE000511.1_2|toxin.fa.cluster196.fa|0.043|16|73|806098|806317|Helicobacter pylori 26695, complete genome   	",
                        "AE000511.1_3|transport.fa.cluster11.fa|2.4e-26|6|189|803756|804368|Helicobacter pylori 26695, complete genome  ",
                        "AE000511.1_3|transport.fa.cluster14.fa|0.0017|17|100|804092|804407|Helicobacter pylori 26695, complete genome  ",
                        "---------- 																									",
                        "AE000512.1_5|transport.fa.cluster4.fa|1.4e-96|18|216|570355|569572|Thermotoga maritima MSB8, complete genome   ",
                        "AE000512.1_5|transport.fa.cluster13.fa|2.5e-67|18|197|570283|569701|Thermotoga maritima MSB8, complete genome  ",
                        "AE000512.1_6|toxin.fa.cluster202.fa|9.9e-07|24|71|597593|597359|Thermotoga maritima MSB8, complete genome  	",
                        "AE000512.1_5|transport.fa.cluster14.fa|5.3e-130|20|103|569947|569635|Thermotoga maritima MSB8, complete genome ",
                        "AE000512.1_5|transport.fa.cluster2.fa|1.2e-54|337|379|570328|570121|Thermotoga maritima MSB8, complete genome  ",
                        "AE000512.1_5|transport.fa.cluster12.fa|2.3e-254|16|218|570295|569671|Thermotoga maritima MSB8, complete genome ",
                        "AE000512.1_5|transport.fa.cluster10.fa|2.4e-274|12|215|570310|569650|Thermotoga maritima MSB8, complete genome ",
                        "AE000512.1_5|transport.fa.cluster3.fa|7.6e-44|285|493|570292|569638|Thermotoga maritima MSB8, complete genome  ",
                        "AE000512.1_5|transport.fa.cluster11.fa|2.1e-234|7|204|570268|569677|Thermotoga maritima MSB8, complete genome  ",
                        "AE000512.1_5|immunity.fa.cluster2.fa|1.3e-159|3|198|570301|569677|Thermotoga maritima MSB8, complete genome    "]

            self.operon_file = "test.operon"
            self.bagel_file = "test.bagel"
            self.out = "test.out"
            open(self.bagel_file,'w').write("\n".join(bagel_entries))
            open(self.operon_file,'w').write("\n".join(operons))
        def tearDown(self):
            os.remove(self.bagel_file)
            os.remove(self.operon_file)
            #os.remove(self.out)
            
        def test1(self):
            bagel,both,detect = bagel_compare_species(self.operon_file,self.bagel_file)
            self.assertEquals(bagel,
                set(["Asticcacaulis excentricus",
                "Streptomyces avermitilis",
                "Streptococcus ratti"]))
            self.assertEquals(both,
                set(["Helicobacter pylori"]))
            self.assertEquals(detect,
                set(["Thermotoga maritima"]))

    unittest.main()

