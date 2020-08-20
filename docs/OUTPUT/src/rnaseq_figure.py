import pymisca.ext as pyext
pd = pyext.pd; np = pyext.np
import os
SRCDIR = os.path.dirname(__file__)

from lazydict import LazyDictionary
rnaseq_figure = template = job = LazyDictionary()

# with pyext.getPathStack([]) as stack:
#     ! ls -lhtr *.pk

from util import _get_file
###############################
#### Data loading #############
rnaseq = rnaseq_raw = pyext.readData(_get_file('/home/feng/envs/Fig_POLYQ/rnaseq.pk'))
# rnaseq = rnaseq_raw = pyext.readData(pyext.f('{SRCDIR}/static.envs.Fig_POLYQ.rnaseq.pk'))
rnaseq = rnaseq.copy()
rnaseq.loc[:] = rnaseq.apply(pyext.log2p1)
job['rnaseq'] = rnaseq


mcurr0 = pyext.readData(_get_file('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv'))
# mcurr0 = pyext.readData( pyext.f('{SRCDIR}/static.results.0318-makeRNA-polyQ.mcurr0.csv') )
mcurr0.columns = mcurr0.columns.str.upper()
mcurr0['DISP_NAME']= pyext.df__format(mcurr0,'{TEMP}-{ZTIME}-{GTYPE}-{index}')
job['datasets_meta'] = mcurr0


keyDF = pyext.readData(pyext.f('{SRCDIR}/key_ath.csv'))
markers = ['LUX']
job['markers_df'] = markers_df = keyDFC = keyDF.query('BioName in %s'%markers)
#####
###############################


# geneDB =pyext.read__buffer(buf,ext='csv')

# markers = ['ATHB2','HFR1','LUX','GI','PRR7']
# markers = ['ATHB2','HFR1','LUX','GI','PRR7']

job['__doc__']=  '''
signature_profile <-- (biological_marker or pca_eigenvector )

signature_profile + signature_dataset --> signature_targets 

signature_targets ~   + visualise_dataset -> [ visualise_boxplot,
                                           visualise_heatmap,
                                           ]
    
signautre_targets + chipseq_targets -> functional_targets
'''

job['order'] = ['ZTIME','GTYPE','TEMP',]

@pyext.setItem(job)
def signature_datasets(self, key, datasets_meta, order):
#     keys = ['ZTIME','GTYPE','TEMP',]
    keys = order
    mcurr0 = datasets_meta
    mcurr = mcurr0
    mcurr = mcurr.query("TEMP in ['22C','27C']")
    mcurr = mcurr.query("RUNID in ['181R','184R','177R']")
    mcurr = mcurr.query("ZTIME in ['ZT12']")
    mcurr = mcurr.query("GTYPE in ['COL','ELF3','ELF3-OX AT']")
    mcurr = mcurr.query("~BNAME.str.contains('adult')")
    mcurr = mcurr.sort_values(keys,)
    # mcurr.index.isin(rnaseq.columns)
    print mcurr.shape
    mcurr[keys]

    signature_datasets = rnaseq.reindex(columns = mcurr.index)
    return signature_datasets


@pyext.setItem(job)
def signature_datasets_meta(self, key, datasets_meta,
                           signature_datasets):
    mcurr = datasets_meta.reindex(signature_datasets.columns)
    return mcurr



@pyext.setItem(job)
def visualise_datasets(self,key, datasets_meta, order):
#     keys = ['ZTIME','GTYPE','TEMP',]
    keys = order

    mcurr0 = datasets_meta
    mcurr = mcurr0
    mcurr = mcurr.query("RUNID in '181R,184R,177R,178R'.split(',')")
    # mcurr = mcurr.query("TEMP in ['22C','27C']")
    mcurr = mcurr.sort_values(keys)

    visualise_datasets = rnaseq.reindex(columns=mcurr.index)
    csv = mcurr0[keys].reindex(visualise_datasets.columns)
    csv["IN_HEATMAP"] = 1
    csv["IN_BOXPLOT"] = 1
    csv.to_csv( "visualise_datasets.csv")
    print csv.shape
    return visualise_datasets


@pyext.setItem(job)
def signature_heatmap(self,key ,signature_datasets, signature_targets):
    vdf = signature_datasets.reindex(signature_targets)
    # fig
    vdf.loc[:] = sutil.meanNorm(vdf)
    vdf = vdf.T
    im  =pyvis.heatmap(vdf,figsize=[12,7])
    return plt.gcf()





@pyext.setItem(job,"plot/signature_profile")
def _func(self,key, 
         datasets_meta,
         signature_datasets,
         
         pyvis,plt,):
    
    mcurr = datasets_meta.reindex(signature_datasets.columns)
#     mcurr['DISP_NAME']= pyext.df__format(mcurr,'{TEMP}-{ZTIME}-{GTYPE}')

    ax = plt.gca()
    vdf = rnaseq.reindex(markers_df.index,columns=mcurr.index)
    vdf.columns = mcurr['DISP_NAME']
    vdf = vdf.T
    # vdf.co
    vdf.plot(xticks=range(len(vdf)),rot='vertical',ax=ax)
    ax.grid(1)
    ax.set_xlabel('dataset')
    ax.set_ylabel('log2(1+TPM_or_CPM)')
    return ax

import synotil.util as sutil 
job['sutil'] = sutil


import pymisca.vis_util as pyvis
plt = pyvis.plt
job['pyvis'] = pyvis
job['plt'] = pyvis.plt
# %matplotlib inline


@pyext.setItem(job)
def signature_profile( self,key,sutil
                      ,signature_datasets,  markers_df):
    vdf = signature_datasets
    vdf = vdf.reindex(index=markers_df.index)
    vdf.loc[:] = sutil.meanNorm(vdf)
    signature_profile = vdf.mean(axis=0)
    
    return signature_profile


job['signature_CUTOFF'] = 0.99
@pyext.setItem(job, 'signature_targets')
def _func(self,key,
         signature_datasets,
         signature_profile,
          signature_CUTOFF,
          pyvis,
#           WORKDIR,
         ):
    vdf = signature_datasets
    signature_score = vdf.dot(signature_profile)
#     ax.set_ylim(0,5000)
    ppf = pyext.dist2ppf(signature_score)
    silent = 1
    if not silent:
        fig,ax = plt.subplots(1,1)
        ax.set_ylabel('signature_score')
        ax.set_xlabel('percentage')
        plt.scatter(ppf,signature_score)
        ax= plt.gca()
        ax.set_xlim(0.95,1.01)
        ax.grid(1)

        pyvis.abline(x0=signature_CUTOFF)
#     CUTOFF = 0.99 ### top 1% of 
    _targets = vdf.index[ppf > signature_CUTOFF]
#     with pyext.getPathStack([WORKDIR,key],force=1) as stack:
    pyext.printlines(_targets,  pyext.f( '_temp-{key}-.it'))
    return _targets


@pyext.setItem(job)
def visualise_heatmap(self,key, visualise_datasets, signature_targets,
                     datasets_meta):
#     pyext.ipd.display(visualise_datasets.head())
#     vdf.columns = mcurr0.reindex(columns)[ keys].reset_index()
#     pyext.ipd.display(vdf.head())

    vdf = visualise_datasets.reindex(signature_targets)
    
    vdf.columns = datasets_meta.reindex(vdf.columns)['DISP_NAME']
    vdf.loc[:] = sutil.meanNorm(vdf)
    vdf = vdf.T
    im = pyvis.df__heatmap(vdf,
                           cname='value',
#                            figsize=[12,15],
    #                        ytick=vdf.index
                      )
    return plt.gcf()

@pyext.setItem(job)
def visualise_boxplot(self,key, 
                      visualise_datasets, signature_targets,
                      datasets_meta,
                      
                      plt,pyvis,sutil,
                     ):
    vdf = visualise_datasets.reindex(signature_targets)
    vdf.columns = datasets_meta.reindex(vdf.columns)['DISP_NAME']
    
    vdf.loc[:] = sutil.meanNorm( vdf )
#     vdf = vdf.T
    fig,ax = plt.subplots(1,1,figsize=[16,7])
    vdf.boxplot(rot='vertical',ax=ax)
    ax.set_ylabel('meanNorm(log2( 1+ CPM_or_TPM))')
#     im = pyvis.df__heatmap(vdf,figsize=[12,15],
# #                        ytick=vdf.index
#                       )

    return fig

@pyext.setItem(job,"export")
def _func(self,key,OUTDIR,
         signature_targets):
    df = pd.DataFrame({"rnaseq_signature_targets":signature_targets})
    with pyext.getPathStack([OUTDIR],force=1):
        df.to_csv("gene_lists_dataframe.csv",index=0)
#     pyext.printlines(signature_targets)
    pass
