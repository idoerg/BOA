""" 
A pipeline to launch parallelize computation on the 
quorum computing cluster
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from collections import defaultdict

import sys
import os,shutil
import site
import argparse
import string
import numpy
import re
import subprocess
import pickle
import tempfile
import glob
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))

import genbank
import blast
import intergene
import genome
import annotation
import bacteriocin
import nbayes
import rforests
import cdhit
import fasta
import hmmer
import quorum
from mafft import MAFFT 
import gff
import interval_filter
import clique_filter
import faa 

class QuorumPipelineHandler(object):
    def __init__(self,
                 rootdir,
                 genome_dir,
                 intergenes,
                 annotated_genes,
                 bacteriocins,
                 bacteriocin_radius,
                 similarity,
                 bac_evalue,
                 training_labels,
                 training_directory,
                 text_database,
                 intermediate,
                 output,
                 numThreads,
                 formatdb,
                 verbose                
                 ):
        
        #Declare global vars
        self.rootdir            =     rootdir   			
        self.genome_dir         =     genome_dir
        self.bacteriocins   	=	  bacteriocins  		
        self.bacteriocin_radius =	  bacteriocin_radius	
        self.similarity	        =	  similarity			
        self.bac_evalue	        =	  bac_evalue
        self.training_labels    =     training_labels			
        self.training_directory =     training_directory
        self.text_database      =     text_database			
        self.intermediate   	=	  intermediate  		
        self.numThreads         =     numThreads
        self.formatdb           =     formatdb				
        self.verbose			=	  verbose
        
        if not os.path.exists(self.intermediate):
            #shutil.rmtree(self.intermediate)
            os.mkdir(self.intermediate)
        if intergenes!=None and annotated_genes!=None: 
            self.intergenes	        =	  intergenes			
            self.annotated_genes    =	  annotated_genes   	
        else:
            self.intergenes = "%s/intergeneDB.fa"%self.intermediate
            self.annotated_genes = "%s/annotated_genesDB.fa"%self.intermediate
        
        self.all_fasta          ="%s/all.fna"%self.intermediate
        self.all_faidx          ="%s/all.fai"%self.intermediate
        self.six_fasta          ="%s/all_trans.fna"%self.intermediate
        self.six_faidx          ="%s/all_trans.fai"%self.intermediate
        self.faa                ="%s/all.faa"%self.intermediate
        self.faaidx             ="%s/all.faaidx"%self.intermediate
        self.gff                ="%s/all.gff"%self.intermediate
        self.blasted_tab_bacteriocins = "%s/blasted_bacteriocins.txt"%self.intermediate
        self.blasted_fasta_bacteriocins = "%s/blasted_bacteriocins.fa"%self.intermediate
        self.cand_context_genes_tab = "%s/cand_context_genes.txt"%self.intermediate
        self.cand_context_genes_fasta = "%s/cand_context_genes.fa"%self.intermediate
        self.cand_context_cluster = "%s/cand_context_cluster"%self.intermediate
        self.blast_context_out = "%s/classify"%self.intermediate
        self.hmmer_context_out = "%s/hmmer"%self.intermediate
        
        #Declare object handlers
        self.clusterer = None                       
        self.textClassifier = None
        self.hmmers = []   
        self.nbpickle = "nb.zip"
        self.clusterpickle = "cluster.zip"
        self.textout = "text_out.txt"
        self.jobs  = []
        self.split_files = [] 
        self.classes = ["toxin","modifier","immunity","transport","regulator"]
        self.class_files = ["%s/%s.fa"%(self.intermediate,s) for s in self.classes]
        self.hmmer_class_out = ["%s/%s.out"%(self.intermediate,s) for s in self.classes] 
        self.operons_out = "%s/operons.txt"%self.intermediate
        self.pred_operons_out = "%s/predicted_operons.txt"%self.intermediate
        
        self.batch_files = []
    def cleanup(self):
        print "Cleaning up"
        for job in self.jobs:
            job.erase_files()
        for fname in self.split_files:
            if os.path.exists(fname): os.remove(fname)
        for hmm in self.hmmers:
            hmm.cleanUp() 
        for i in xrange(len(self.class_files)):
            if os.path.exists(self.class_files[i]): os.remove(self.class_files[i])
        print glob.glob("tmp*")
        for file in glob.glob("tmp*"):
            os.remove(file)
        for file in self.batch_files:
            for extra in glob.glob("%s*"%file):
                os.remove(extra)
            
    """ Builds database such as the 
        intergenic database
        annotated genes database 
        naive bayes model"""
    def preprocess(self):
        print "Preprocessing"
        #Combine all genome files into a single genome fasta file
        fasta.go(self.genome_dir,
                 self.all_fasta,
                 self.all_faidx,
                 self.six_fasta,
                 self.six_faidx) 
        indexer = fasta.Indexer(self.all_fasta,self.all_faidx)
        indexer.index()
        indexer.load()
        intergene.go(self.genome_dir,self.intergenes)
        annotation.go(self.genome_dir,self.annotated_genes,index_obj=indexer) 
        #Combine all gff files together
        outhandle = open(self.gff,'w')
        for root, subFolders, files in os.walk(self.genome_dir):
            for fname in files:
                genome_files = []
                organism,ext = os.path.splitext(os.path.basename(fname))
                absfile=os.path.join(root,fname)
                if ext==".gff":
                    shutil.copyfileobj(open(absfile),outhandle)
        outhandle.close()
        
        tmpfile = "tmp%d.faa"%(os.getpid())
        outhandle = open(tmpfile,'w')
        for root, subFolders, files in os.walk(self.genome_dir):
            for fname in files:
                genome_files = []
                organism,ext = os.path.splitext(os.path.basename(fname))
                absfile=os.path.join(root,fname)
                if ext==".faa":
                    shutil.copyfileobj(open(absfile),outhandle)
        outhandle.close()
        faa.reformat(tmpfile,self.faa)
        os.remove(tmpfile)
        
        faaindex = fasta.Indexer(self.faa,self.faaidx)
        faaindex.index()
        
    """  Runs blast to identify bacteriocins and context genes"""
    def blast(self,njobs=1):
        print "Blasting"
        
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_bacfiles = ["%s/bacteriocin.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        self.split_files += split_bacfiles
        split_bachandles = [open(f,'w') for f in split_bacfiles]
        out_fnames = ["%s/blasted.%d"%(self.intermediate,i) for i in xrange(njobs)]
        out_bac = ["%s.bacteriocins.txt"%(out) for out in out_fnames]
        out_genes = ["%s.annotated.txt"%(out) for out in out_fnames]
        index=0
        for record in SeqIO.parse(self.bacteriocins,"fasta"):
            split_bachandles[index].write(">%s\n%s\n"%(str(record.id),
                                                       str(record.seq)))
            index=(index+1)%njobs
        #Close files
        for handle in split_bachandles: handle.close()
        print "evalue",self.bac_evalue
        if self.formatdb: 
            blast_cmd = ' '.join([
                        """module load anaconda; module load blast;module load blast+;"""
                        """python %s/src/genome/bacteriocin.py  """,
                        """ --genome-files %s        """,
                        """ --annotated-genes=%s     """,
                        """ --intergenes=%s          """,
                        """ --bacteriocins=%s   	 """,
                        """ --bacteriocin-radius=%d  """,
                        """ --bac-evalue=%s 		 """,
                        """ --num-threads=%d    	 """,
                        """ --intermediate=%s   	 """,
                        """ --output=%s     		 """,
                        """ --formatdb               """,
                        """ --verbose                """
                        ])
        else:
            blast_cmd = ' '.join([
                        """module load anaconda; module load blast;module load blast+;"""
                        """python %s/src/genome/bacteriocin.py  """,
                        """ --genome-files %s        """,
                        """ --annotated-genes=%s     """,
                        """ --intergenes=%s          """,
                        """ --bacteriocins=%s   	 """,
                        """ --bacteriocin-radius=%d  """,
                        """ --bac-evalue=%s 		 """,
                        """ --num-threads=%d    	 """,
                        """ --intermediate=%s   	 """,
                        """ --output=%s     		 """,
                        """ --verbose                """
                        ])
            
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = blast_cmd%(self.rootdir,
                            self.all_fasta,
                            self.annotated_genes,
                            self.intergenes,
                            split_bacfiles[i],
                            self.bacteriocin_radius,
                            str(self.bac_evalue),
                            self.numThreads,
                            self.intermediate,
                            out_fnames[i])
            
            batch_file = "%s/blast%i.%d.job"%(os.getcwd(),i,os.getpid())
            self.batch_files.append(batch_file)
            proc = quorum.Popen(cmd,shell=True,batch_file=batch_file,
                                stdin=quorum.PIPE,stdout=quorum.PIPE ,threads=self.numThreads)
            proc.submit()
            #proc.output = out_fnames[i]
            jobs.append(proc)
            self.jobs.append(proc)
        
        for job in jobs: job.wait() 
        
        """ Collect all of the results from the jobs"""
        bacteriocins_out = open(self.blasted_fasta_bacteriocins,'w')
        context_genes_out = open(self.cand_context_genes_fasta,'w')
        out_bac = ["%s.bacteriocins.txt"%(out) for out in out_fnames]
        out_genes = ["%s.annotated.txt"%(out) for out in out_fnames]
        for i in xrange(njobs):
            if os.path.exists(out_bac[i]):
                shutil.copyfileobj(open(out_bac[i]),bacteriocins_out)
            if os.path.exists(out_genes[i]):
                shutil.copyfileobj(open(out_genes[i]),context_genes_out)
        bacteriocins_out.close()
        context_genes_out.close()
        
        
    """ Clusters bacteriocins and context genes together"""
    def cluster(self,preprocess=False,numThreads=8,mem=6000):
        print "Clustering"
        if preprocess:
            fasta.preprocess(self.cand_context_genes_tab,
                            self.cand_context_genes_fasta)
            
        
        self.clusterer = cdhit.CDHit(self.cand_context_genes_fasta,
                                     self.cand_context_cluster,
                                     self.similarity,
                                     numThreads,mem)
        self.clusterer.run()
        self.clusterer.parseClusters()
        self.clusterer.dump(self.clusterpickle)
    """ Identifies context genes using BLAST"""
    def blastContextGenes(self,njobs=2):
        print "Blasting Context Genes"
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_fastafiles = ["%s/context.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        #self.split_files += split_fastafiles
        split_fastahandles = [open(f,'w') for f in split_fastafiles]
        out_classes = ["%s/contextout.%d"%(self.intermediate,i) for i in xrange(njobs)]
        index=0
        for record in SeqIO.parse(self.cand_context_genes_fasta,"fasta"):
            if len(record.seq)<=1: continue #To weed out weird entries
            split_fastahandles[index].write(">%s\n%s\n"%(str(record.id),
                                                         fasta.format(str(record.seq))))
            index=(index+1)%njobs
        #Close files
        for handle in split_fastahandles: handle.close()
        context_cmd = ' '.join([
                                 """module load anaconda; module load blast;module load blast+;""",
                                 """python %s/src/genome/context_gene.py""",
                                 """--training-directory=%s""",
                                 """--training-labels=%s""",
                                 """--query=%s""",
                                 """--intermediate=%s""",
                                 """--num-threads=%d""",
                                 """--output=%s"""         
                                 ])    
        
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = context_cmd%(self.rootdir,
                               self.training_directory,
                               self.training_labels,
                               split_fastafiles[i],
                               self.intermediate,
                               self.numThreads,
                               out_classes[i]
                               )
            
            batch_file = "%s/context_blast%i.%d.job"%(os.getcwd(),i,os.getpid())
            self.batch_files.append(batch_file)
            
            proc = quorum.Popen( cmd,shell=True,batch_file=batch_file,
                                 stdin=quorum.PIPE,stdout=quorum.PIPE ,threads=self.numThreads) 
            proc.submit()
            #proc.output = out_classes[i]
            jobs.append(proc)
            self.jobs.append(proc)
        for job in jobs: job.wait()
        """ Collect all of the results from the jobs"""
        context_out = open(self.blast_context_out,'w')
        for i in xrange(njobs):
            if os.path.exists(out_classes[i]):
                shutil.copyfileobj(open(out_classes[i]),context_out)
        context_out.close()    
        pass
    
    """ 
    Creates HMMER profiles of all of the blasted context genes and the 
    bacteriocins to find more candidate context genes and bacteriocins
    """
    def hmmerGenes(self,min_cluster=10,njobs=2):
        print "Hmmering Genes"
        # Aggregate all of the bacteriocins and context genes
        # into 5 files representing each class 
        class_handles = [open(c,'w') for c in self.class_files ]
        for record in SeqIO.parse(self.blast_context_out,"fasta"):
            label = record.id.split('|')[-1]
            index = self.classes.index(label)
            SeqIO.write(record,class_handles[index],"fasta")
        #Move all bacteriocins into toxins file
        shutil.copyfileobj(open(self.blasted_fasta_bacteriocins),
                           class_handles[0])
        for c in class_handles: c.close()
        for fname in self.class_files:
            assert os.path.getsize(fname)>0
        #Remove all duplicate ids from this files
        for fname in self.class_files:
            f = "%s.%d"%(fname,os.getpid())
            fasta.remove_duplicates(fname,f)
            os.rename(f,fname)
        for fname in self.class_files:
            assert os.path.getsize(fname)>0
        #self.class_files = self.class_files[::-1] #Reverse list for testing]
        # Run separate instances of HMMER for each class
        print "threads",self.numThreads
        self.hmmers = [hmmer.HMMER(f,quorum,self.numThreads,min_cluster) 
                       for f in self.class_files]
        for i in xrange(len(self.class_files)):
            H = self.hmmers[i]
            H.writeClusters(similarity=0.7,memory=3000)
            H.HMMspawn(msa=MAFFT,njobs=njobs)
            H.search(self.six_fasta,self.hmmer_class_out[i],njobs)
    
    """ Writes clusters and their corresponding sequences"""
    def writeClusters(self,clusters,gff_index,faaindex,outhandle):
        clusternum=0
        for cluster in clusters:
            genes = []
            for gene in cluster:
                acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=gene.split("|")
                genes.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description))
            orfs = gff_index.translate_orfs(genes,faaindex)
            for orf in orfs:
                id,seq = orf
                acc,clrname,full_evalue,env_st,env_end,protid = id
                function = clrname.split('.')[0]
                outhandle.write(">accession=%s|function=%s|start=%s|end=%s|evalue=%s|protid=%s|cluster_%d\n%s\n"%
                                (acc,function,env_st,env_end,str(full_evalue),protid,clusternum,
                                 fasta.format(seq))) 
            clusternum+=1

    """ Finds operons by constructing graphs and finding cliques 
    TODO: Move these parameters to main pipeline handler object"""
    def cliqueFilter(self,clique_radius=50000,functions = ["toxin","modifier","immunity","transport","regulator"]):
        print "Clique filtering","Looking for cliques with",functions
        
        toxin_hits     = hmmer.parse("%s/toxin.out"%self.intermediate)
        modifier_hits  = hmmer.parse("%s/modifier.out"%self.intermediate)
        immunity_hits  = hmmer.parse("%s/immunity.out"%self.intermediate)
        regulator_hits = hmmer.parse("%s/regulator.out"%self.intermediate)
        transport_hits = hmmer.parse("%s/transport.out"%self.intermediate)
                
        genefile = gff.GFF(gff_file=self.gff,fasta_index=self.six_faidx)
        toxin_hits     = genefile.call_orfs(toxin_hits    )
        modifier_hits  = genefile.call_orfs(modifier_hits )
        immunity_hits  = genefile.call_orfs(immunity_hits )
        regulator_hits = genefile.call_orfs(regulator_hits)
        transport_hits = genefile.call_orfs(transport_hits)
        all_hits = toxin_hits+modifier_hits+immunity_hits+regulator_hits+transport_hits
        #Find operons with at least a toxin and a transport
        all_hits = interval_filter.unique(all_hits)
        # #Sort by start/end position and genome name
        all_hits=sorted(all_hits,key=lambda x: x[6])   
        all_hits=sorted(all_hits,key=lambda x: x[5])
        all_hits=sorted(all_hits,key=lambda x: x[0])
        all_hits=sorted(all_hits,key=lambda x: x[-1])
        print "Size of All hits",sys.getsizeof(all_hits)
        print "Size of genefile",sys.getsizeof(genefile)
        del toxin_hits
        del modifier_hits
        del immunity_hits
        del regulator_hits
        del transport_hits
        
        clusters = clique_filter.findContextGeneClusters(all_hits,self.six_faidx,
                                                         radius=clique_radius,
                                                         backtrans=False,
                                                         functions=["toxin","transport"])
         
        faaindex = fasta.Indexer(self.faa,self.faaidx)
        faaindex.index()
        faaindex.load()
        
        print "Size of FAA",sys.getsizeof(faaindex)
        
        outhandle = open(self.operons_out,'w')
        self.writeClusters(clusters,genefile,faaindex,outhandle) 
        outhandle.close()
        #Predict operons based on just context genes
        clusters = clique_filter.findContextGeneClusters(all_hits,self.six_faidx,
                                                         radius=clique_radius,
                                                         backtrans=False,   
                                                         functions=["modifier","regulator","immunity","transport"])
        
        outhandle = open(self.pred_operons_out,'w')
        self.writeClusters(clusters,genefile,faaindex,outhandle) 
        outhandle.close()
            
        #=======================================================================
        # toxin_hits     = parse("%s/toxin.out"%self.intermediate)
        # modifier_hits  = parse("%s/modifier.out"%self.intermediate)
        # immunity_hits  = parse("%s/immunity.out"%self.intermediate)
        # regulator_hits = parse("%s/regulator.out"%self.intermediate)
        # transport_hits = parse("%s/transport.out"%self.intermediate)
        # all_hits = toxin_hits+modifier_hits+immunity_hits+regulator_hits+transport_hits 
        # 
        # #Sort by start/end position
        # all_hits=sorted(all_hits,key=lambda x: x[6])   
        # all_hits=sorted(all_hits,key=lambda x: x[5])
        # 
        # #Sort by genome name
        # all_hits=sorted(all_hits,key=lambda x: x[0])
        # clique_cmd = ' '.join([
        #                 """module load anaconda; module load blast;module load blast+;"""
        #                 """python %s/src/filter/clique_filter.py  """,
        #                 """ --input=%s               """,
        #                 """ --functions %s           """,
        #                 """ --clique-radius=%d       """,
        #                 """ --faidx=%s               """]) 
        # prevGenome = None
        # jobs,buf,clusters = [],[],[]
        # i = 0
        # for hit in all_hits:
        #     if prevGenome == None:      
        #         prevGenome = hit[-1]
        #     elif prevGenome == hit[-1]: 
        #         buf.append(hit)
        #     else: 
        #         #Submit job
        #         bufstr = hmmerstr(buf)
        #         
        #         split_file = "%s/hmm.%d.%d.clique.txt"%(self.intermediate,i,os.getpid())
        #         with open(split_file,'w') as outhandle: outhandle.write(bufstr)
        #         batch_file = "%s/clique.%d.%d.job"%(os.getcwd(),i,os.getpid())
        #         self.split_files.append(split_file)
        #         self.batch_files.append(batch_file)
        #         cmd = clique_cmd%(self.rootdir,
        #                           split_file,
        #                           ' '.join(functions),
        #                           clique_radius,
        #                           self.faidx)
        #         proc = quorum.Popen(cmd,shell=True,batch_file=batch_file,
        #                             stdin=quorum.PIPE,stdout=quorum.PIPE)
        #         proc.submit()
        #         jobs.append(proc)
        #         self.jobs.append(proc)
        #         buf = [hit]
        #         prevGenome = hit[-1]    
        #         i+=1
        #         
        #         
        # operon_strings = []
        # for job in jobs: 
        #     job.wait()
        #     stdout_str = job.ofile_string()
        #     operon_strings.append(stdout_str)
        # 
        # open(self.operons_out,'w').write(''.join(operon_strings))
        #=======================================================================
        
    """ Deprecated: Classifies individual bacteriocins and context genes based on their text"""
    def textmine(self,njobs=1,dbtype="sqlite3",numTrees=1000):
        index=0
        print "Classifying"
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_fastafiles = ["%s/cluster.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        self.split_files += split_fastafiles
        split_fastahandles = [open(f,'w') for f in split_fastafiles]
        out_classes = ["%s/classify.%d"%(self.intermediate,i) for i in xrange(njobs)]
        for record in SeqIO.parse(self.cand_context_cluster,"fasta"):
            split_fastahandles[index].write(">%s\n%s\n"%(str(record.id),
                                                       str(record.seq)))
            index=(index+1)%njobs
        #Close files
        for handle in split_fastahandles: handle.close()
        classify_cmd = ' '.join([
                                 """module load anaconda; module load blast;""",
                                 """python %s/src/classify/rforests.py""",
                                 """--training-directory=%s""",
                                 """--training-labels=%s""",
                                 """--text-database=%s""",
                                 """--database-type=%s""",
                                 """--num-trees=%d""",
                                 """--fasta=%s""",
                                 """--output=%s"""         
                                 ])    
        
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = classify_cmd%(self.rootdir,
                                self.training_directory,
                                self.training_labels,
                                self.text_database,
                                dbtype,
                                numTrees,
                                split_fastafiles[i],
                                out_classes[i]
                                )
            
            batch_file = "%s/classify%i.%d.job"%(os.getcwd(),i,os.getpid())
            self.batch_files.append(batch_file)
            proc = quorum.Popen( cmd,shell=True,batch_file=batch_file,
                                 stdin=quorum.PIPE,stdout=quorum.PIPE ,threads=self.numThreads) 
            proc.submit()
            #proc.output = out_classes[i]
            jobs.append(proc)
            self.jobs.append(proc)
        
        for job in jobs: job.wait()
        
        """ Collect all of the results from the jobs"""
        classify_out = open(self.blast_context_out,'w')
        for i in xrange(njobs):
            shutil.copyfileobj(open(out_classes[i]),classify_out)
        classify_out.close()
        
        
        
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins and context genes')
    parser.add_argument(\
        '--pipeline-section', type=str, required=False, default="all",
        help='Section of the pipeline to run (all, preprocess, blast, context, hmmer, clique)')
    parser.add_argument(\
        '--root-dir',type=str, required=False,
        help='Root directory')
    parser.add_argument(\
        '--genome-dir', type=str, required=False,
        help='Directory of all of the fasta/genbank files')
    parser.add_argument(\
        '--intergenes', type=str, required=False,default=None,
        help='FASTA files containing intergenic regions (default=None)')
    parser.add_argument(\
        '--annotated-genes', type=str, required=False,default=None,
        help='FASTA files containing annotated genetic regions (default=None)')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,default=None,
        help='The bacteriocin proteins that are to be blasted')
    parser.add_argument(\
        '--bacteriocin-radius', type=int, required=False, default=50000,
        help='The search radius around every specified bacteriocin')
    parser.add_argument(\
        '--similarity', type=float, required=False, default=0.7,
        help='Clustering similarity')    
    parser.add_argument(\
        '--cluster-size', type=int, required=False, default=10,
        help='Filters all clusters below this threshold')    
    parser.add_argument(\
        '--bac-evalue', type=float, required=False, default=0.00001,
        help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--functions', type=str, nargs="+", default=None, required=False,
        help='The list of functions to look for (e.g. toxin, modifier, transport, immunity, regulator)')
    parser.add_argument(\
        '--training-labels',type=str,required=False,default=None,
        help='Training labels for Naive Bayes')
    parser.add_argument(\
        '--training-directory',type=str,required=False,default=None,
        help='''Training directory containing all 
                of the context genes required for finding context genes''')
    parser.add_argument(\
        '--text-database', type=str, required=False,
        help='SQL database containing text annotations (default=None)')
    parser.add_argument(\
        '--intermediate', type=str, required=False,default='.',
        help='Directory for storing intermediate files')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file basename for filtered annotationed regions and bacteriocins')
    parser.add_argument(\
        '--formatdb', action='store_const', const=True, default=False,
        help='Indicates if formatdb should be run')
    parser.add_argument(\
        '--num-threads', type=int, required=False, default=1,
        help='The number of threads to be run by BLAST')
    parser.add_argument(\
        '--num-jobs', type=int, required=False, default=1,
        help='The number of parallel jobs to be run')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Messages for debugging')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    print "Intermediate",args.intermediate
    print "evalue",args.bac_evalue
    if not args.test:
        if args.functions==None:
            functions = ["toxin","modifier","immunity","transport","regulator"]
        proc = QuorumPipelineHandler(args.root_dir,
                                     args.genome_dir,
                                     args.intergenes,
                                     args.annotated_genes,
                                     args.bacteriocins,
                                     args.bacteriocin_radius,
                                     args.similarity,
                                     args.bac_evalue,
                                     args.training_labels,
                                     args.training_directory,
                                     args.text_database,
                                     args.intermediate,
                                     args.output,
                                     args.num_threads,
                                     args.formatdb,
                                     args.verbose                
                                    )                                     
        if args.pipeline_section=="all":
            proc.preprocess()
            proc.blast(njobs=args.num_jobs)
            proc.blastContextGenes(njobs=args.num_jobs)
            proc.hmmerGenes(args.cluster_size,args.num_jobs)
            proc.cliqueFilter(args.bacteriocin_radius,functions=args.functions)            
        elif args.pipeline_section=="blast":
            proc.blast(njobs=args.num_jobs)
            proc.blastContextGenes(njobs=args.num_jobs)
            proc.hmmerGenes(args.cluster_size,args.num_jobs)
            proc.cliqueFilter(args.bacteriocin_radius,functions=args.functions)
        elif args.pipeline_section=="context":
            proc.blastContextGenes(njobs=args.num_jobs)
            proc.hmmerGenes(args.cluster_size,args.num_jobs)
            proc.cliqueFilter(args.bacteriocin_radius,functions=args.functions)
        elif args.pipeline_section=="hmmer":
            proc.hmmerGenes(args.cluster_size,args.num_jobs)
            proc.cliqueFilter(args.bacteriocin_radius,functions=args.functions)
        elif args.pipeline_section=="clique":
            proc.cliqueFilter(args.bacteriocin_radius,args.functions)
            pass
        #from time import sleep
        #sleep(100)
        #proc.cleanup()
    else:
        del sys.argv[1:]
        import unittest
        import test_modules
        class TestRun(unittest.TestCase):
            def setUp(self):
                self.root = os.environ['BACFINDER_HOME']
                self.exampledir = "%s/example/Streptococcus_pyogenes"%self.root
                self.bacdir = "%s/bacteriocins"%self.root
                self.genome_files = test_modules.getFNA(self.exampledir) 
                self.six_frame_genome = "%s/example/all_orfs.fna"%self.root
                self.six_frame_genome_index = "%s/example/all_orfs.fai"%self.root
                
                self.bacteriocins = "%s/bagel.fa"%self.bacdir
                self.intergenes = "test_intergenes.fa"
                self.annotated_genes = "test_genes.fa"
                self.textdb = "%s/db/bacteria_database"%(self.root)
                self.textpickle = "%s/db/bacteria_table.zip"%(self.root)
                
                #self.intergenes = "%s/db/intergenes.txt"%(self.root)
                #self.annotated_genes = "%s/db/annotated_genes.txt"%(self.root)
                self.intermediate = "intermediate"
                #self.training_labels = "%s/data/training/training.txt"%self.root
                self.training_directory = "%s/data/training/protein"%self.root
                self.training_labels = "%s/data/training/training_proteins.txt"%self.root
                if not os.path.exists(self.intermediate):
                    os.mkdir(self.intermediate)
                self.bac_evalue = 0.000001
                self.formatdb = True
                self.bacteriocin_radius = 50000
                self.verbose = True
                self.keep_tmp = False
                self.similarity = 0.65
                self.numThreads = 6
                self.output = "out"
                self.keep_tmp = True
                self.proc = None
            def tearDown(self):
                #self.proc.cleanup()
                pass
            def testrun(self):
                print "Test Run"
                self.proc = QuorumPipelineHandler(                       
                                         self.root,
                                         self.exampledir,
                                         self.intergenes,
                                         self.annotated_genes,
                                         self.bacteriocins,
                                         self.bacteriocin_radius,
                                         self.similarity,
                                         self.bac_evalue,
                                         self.training_labels,
                                         self.training_directory,
                                         self.textdb,
                                         self.intermediate,
                                         self.output,
                                         self.numThreads,
                                         self.formatdb,
                                         self.verbose                
                                        ) 
                
                self.proc.preprocess()
                self.assertTrue(os.path.getsize(self.annotated_genes) > 0)
                self.assertTrue(os.path.getsize(self.intergenes) > 0)
                self.assertTrue(os.path.getsize(self.proc.all_fasta) > 0)
                self.assertTrue(os.path.getsize(self.proc.six_fasta) > 0)
                self.assertTrue(os.path.getsize(self.proc.six_faidx) > 0)
                
                self.proc.blast(njobs=10)
                
                self.assertTrue(os.path.getsize(self.proc.blasted_fasta_bacteriocins) > 0)
                self.assertTrue(os.path.getsize(self.proc.cand_context_genes_fasta) > 0)
                
                self.proc.blastContextGenes(njobs=10)
                self.assertTrue(os.path.getsize( self.proc.blast_context_out ) > 0)
                
                self.proc.hmmerGenes(min_cluster=1,njobs=8)
                for fname in self.proc.class_files:
                    self.assertTrue(os.path.getsize(fname)>0)
                
                self.proc.cliqueFilter(clique_radius=100000)
                self.assertTrue(os.path.getsize(self.proc.operons_out)>0)
                self.assertTrue(os.path.getsize(self.proc.pred_operons_out)>0)
        unittest.main()       
        
        
        
        
        
