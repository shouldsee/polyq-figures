import pymisca.ext as pyext
pd = pyext.pd; np = pyext.np

import pymisca.vis_util as pyvis
plt = pyvis.plt
import path,sys,os

import synotil.util as sutil
import pymisca.util as pyutil

from lazydict import LazyDictionary
template = job = LazyDictionary()
__FILE__ = path.Path(__file__)
DEPENDS = {
'meta_chip':
    (pyext.readData, '/home/feng/meta/meta_chip.tsv'),
'fig_meta':
    (pyext.readData, __FILE__.dirname()/'0726-figure-meta.tsv'),
'chipseq_targets_peaks_file':
    ( None, '/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak'),
    
}

for k in DEPENDS:
    v = DEPENDS[k]
    assert    pyext.file__notEmpty(v[1]),(v[1],)
    if v[0] is not None:
        job[k] = v[0](v[1]) 
    else:
        job[k] = v[1]
    

@pyext.setItem(job, "chipseq_targets_peaks")
def _func(self,key,):
    
    pass

@pyext.setItem(job, "chipseq_targets_genes")
def _func(self, key, chipseq_targets_peaks):
    pass


@pyext.setItem(job, "chipseq_pileup_boxplot")
def _func(self,key, 
          OUTDIR,
          fig_meta, 
          meta_chip,
          chipseq_targets_peaks_file,
         ):
    from pymisca.plotters import plotters
    plotters.fig__save    
#     key = 'figS4E_o'
    key = 'figS4E_0905'
    df = fig_meta[[ key ]].dropna().astype(int)
    df = df.sort_values( key )
    DATA_ACC_LIST = df.index.tolist()
    mcurr = meta_chip.reindex( DATA_ACC_LIST )
    bwFiles = mcurr['RPKMFile']
    
    # bedFile = '/home/feng/envs/Fig_POLYQ/bedFile.bed'
#     bedFile = '/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak'
    res = sutil.extract_bigwig_multiple(
        bwFiles=mcurr.RPKMFile,
        bedFile=chipseq_targets_peaks_file,
#         radius=100,stepSize=10,NCORE=6,
        radius=100,stepSize=10,NCORE=6,
        outIndex=pyext.df__format(mcurr,'{bname}_{index}'),)
    
    tab = dfc = pyutil.colGroupMean(res,).apply(sutil.log2p1)

    with pyext.getPathStack([OUTDIR],force=1):
        plt.figure(figsize=[12,8])
        dfc.boxplot(rot='vertical',)
        plt.ylabel('average binding in RPKM')
#         plt.gcf().savefig('fig1.png')
#         plt.gcf().
        plotters.fig__save(plt.gcf(), pyext.f('{key}.png'))


        # res.head()
        tab2 = res.mean(axis=0)
        pile = tab2.to_frame().reset_index().pivot_table(index='bwFile',columns='pos',values=0)
        pile = pile.T
        pyvis.df__heatmap(pile, figsize=[12,12])
#         pyvis.heatmap(pile, xtick=pile.columns,ytick=pile.index)
        plotters.fig__save(plt.gcf(), pyext.f('{key}_pileup_heatmap.png') )
#         plt.gcf().savefig(pyext.f('{key}_pileup_heatmap.png'))

import matplotlib.pyplot as plt
import synotil.qcplots

import matplotlib
# matplotlib.style.use("ggplot")
# %matplotlib inline
@pyext.setItem(job, "chipseq_targets_genes_job")
def _func(job,key,WORKDIR,chipseq_targets_peaks):
    fig,axs = plt.subplots(2,2,figsize=[12,12])
    axs= axs.ravel()

    with pyext.getPathStack([WORKDIR,key],force=1):
        GSIZE = "/home/feng/ref/ATH-TAIR10/genome.sizes"
        peakFile = chipseq_targets_peaks['LAST_FILE']
        featFile = sutil.bed__leftSummit("/home/feng/ref/ATH-TAIR10/annotation/genes.gtf.cds",GSIZE=GSIZE)
        pyext.file__link(peakFile,"PEAK_FILE.bed",force=1)
        pyext.file__link(featFile,"FEAT_FILE.bed",force=1)
        res = synotil.qcplots.qc_summitDist(
#             chipseq_targets_genes_peaks,
            peakFile,
            featFile,
            GSIZE=GSIZE,
            axs=axs,
            CUTOFF=500,
        )
        df = res[0]
#         pyext.readData(peakFile,header=None).set_index(3,drop=0).reindex(df['acc'])
#         df =  df.set_index('acc',drop=0).reindex(pyext.readData(peakFile,header=None).set_index(3))
        df.to_csv("PEAK_DIST.csv",index=0)

        lst = df.dropna()['feat_acc'].drop_duplicates().values.tolist()
        with pyext.getPathStack([WORKDIR,key],force=1) as stack:
            pyext.printlines(lst,"OUT.it")        
            
    return dict(LAST_DIR = WORKDIR / key)

@pyext.setItem(job,"chipseq_targets_peaks")
def _func(job,key,WORKDIR):
    GSIZE = "/home/feng/ref/ATH-TAIR10/genome.sizes"
    with pyext.getPathStack([WORKDIR,key],force=1) as stack:
        res = sutil.bed__summit('/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',GSIZE=GSIZE,inplace=False)
        res = stack.d/res
    return dict(LAST_FILE=res)

@pyext.setItem(job, "chipseq_targets_genes")
def _func(job,key, chipseq_targets_genes_job):
    return list(pyext.readData(
    chipseq_targets_genes_job['LAST_DIR']/"OUT.it"
    ))


@pyext.setItem(job,"export")
def _func(self,key,
          chipseq_pileup_boxplot,
          chipseq_targets_genes,
         chipseq_targets_peaks,):
    pass
