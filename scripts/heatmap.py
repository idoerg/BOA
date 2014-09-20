###./imp.py deletedOrg /home/asmariyaz/Desktop/phylo_order /home/asmariyaz/Desktop/txtnames
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cbook as cbook
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage 
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from Bio import Phylo




def produce_tree():
    tree=Phylo.read("/home/asmariyaz/Desktop/reorder.nwk","newick")
    Phylo.draw(tree, axes=phyl_ax)
    tree_f=plt.gcf()
    print type(tree_f)
    return tree_f

""" Lookup species names in tree """
def lookup_by_names(tree):
    names = {}
    for clade in tree.get_terminals():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names
        
""" Remove all leaves that aren't contained in names """
def prune_tree(tree,target_species):
    leafs = tree.get_terminals()
    all_species = set([l.code for l in leafs])
    target_species = set(target_species)
    outlier_species = all_species.difference(target_species)
    name_dict = lookup_by_names(tree)
    for outlier in outlier_species:
        tree.prune(name_dict[outlier])
    return tree

def produce_hm(full_len,outlier,max_num,min_num,file_handle,txtnames):
    tree=Phylo.read("/home/asmariyaz/Desktop/reorder.nwk","newick")
    a=1.0
    data=np.array(full_len)
    cmap = mpl.cm.hot
    if outlier==0:
       threshold=1
    else:
       threshold=outlier-0.01
    
    fig = plt.figure(figsize=(25,25))
    plt.suptitle(file_handle.replace('_del.csv',''),fontsize=22)
    cmap.set_over('green')
    cmap.set_under('grey')
    gs=gridspec.GridSpec(1, 2,height_ratios=[1,1,-2,2] ,width_ratios=[1,1,-2,2],hspace=0,wspace=0) 
    phyl_ax=plt.subplot(gs[0])
    Phylo.draw(tree, axes=phyl_ax,do_show=False)    
    ht_ax=plt.subplot(gs[1])
    ht_ax.set_xlim(0,34)
    ht_ax.set_ylim(0,34)
    ht_ax.grid(True, which='both')
    
    ##cb_ax,kw =mpl.colorbar.make_axes(ht_ax, shrink=0.65)
    
    plt.setp(phyl_ax.get_xticklabels(),visible=False)
    plt.setp(phyl_ax.get_yticklabels(),visible=False)
    plt.setp(ht_ax.get_xticklabels(),visible=True)
    plt.setp(ht_ax.get_yticklabels(),visible=False)
    plt.setp(phyl_ax.get_xticklines(),visible=False)
    plt.setp(phyl_ax.get_yticklines(),visible=False)
    plt.setp(ht_ax.get_xticklines(),visible=False)
    plt.setp(ht_ax.get_yticklines(),visible=False)
    xticks=range(34)
    ht_ax.xaxis.set_ticks(xticks)
    ht_ax.yaxis.set_ticks(xticks)
    ht_ax.set_xticklabels(txtnames,rotation=45,fontsize=12.5,alpha=a)
    ht_ax.xaxis.set_tick_params(pad=4)

    for tick in ht_ax.xaxis.get_major_ticks():
        tick.label1.set_horizontalalignment('right')
    
    ##ticks_at = [-2,0,2]
    img = ht_ax.imshow(data, cmap=cmap, interpolation='none',vmax=threshold,aspect='auto',extent=[34,0,34,0],origin='lower')
    
    ##cb = mpl.colorbar.ColorbarBase(ax=cb_ax,cmap=cmap,ticks=ticks_at,extend='neither',**kw)
    ##cb.cmap.set_over('green')

    ##img= mpimg.imread('/home/asmariyaz/Desktop/mytree.png')
    ##phyl_ax.imshow(img,interpolation='bilinear',aspect='auto')
    heatmap_file=fig.savefig('/home/asmariyaz/Desktop/heatmap1/'+file_handle.replace('.csv','')+'.pdf',bbox_inches='tight',dpi=150)
    return heatmap_file



    
##def produce_hm(org_name_list,full_len,outlier,max_num,min_num,file_handle):
##    d=full_len
##    data=np.array(d)
##    fig,ax = plt.subplots()
##    cmap =plt.cm.get_cmap('hot')
##    negative_outlier=0-outlier
##    norm = mpl.colors.Normalize(vmin=negative_outlier,vmax=outlier)
##    
##  
##    plt.xlim(0,35)
##    plt.ylim(0,35)
##    img=plt.imshow(data, interpolation='none',norm=norm, cmap=cmap,vmax=outlier)
##    
##    cb_ax=fig.add_axes([0.85, 0.1, 0.03, 0.8])
##    
##    cb=mpl.colorbar.ColorbarBase(cb_ax,cmap=cmap,norm=norm,extend='both',spacing='uniform')
##    cmap.set_over('green')
##    cmap.set_under('green')
##
##    heatmap_file=fig.savefig('/home/asmariyaz/Desktop/heatmap/'+file_handle+'.png')
##    return heatmap_file




    
##mpl.colorbar.ColorbarBase(ax1,cmap=cmap,orientation='vertical')
  ##cmap.set_bad('green')
##img=plt.imshow(np.ma.masked_values(data, outlier),norm=norm, interpolation='none', cmap=cmap)
