import os,sys
import pymisca.ext as pyext
pd = pyext.pd; np = pyext.np
from util import _get_file
from util import Inputer

# from lazydict import LazyDictionary
# rnaseq_figure = template = job = LazyDictionary()

def fig__saveFile(fig, ofname, transparent = False, **kwargs):
    fig.savefig(ofname,
                bbox_inches='tight',
                        transparent= transparent,
                        facecolor=fig.get_facecolor(),
                        **kwargs)
    return 1


class base(object):
    def __init__(self):
        self.to_override = {}
        # self.input_files = {}
    
    def is_ready(self):
        return not len(self.to_override)




class rnaseq_figure(base):
    __doc__ = '''
    signature_profile <-- (biological_marker or pca_eigenvector )

    signature_profile + signature_dataset --> signature_targets 

    signature_targets ~   + visualise_dataset -> [ visualise_boxplot,
                                               visualise_heatmap,
                                               ]
        
    signautre_targets + chipseq_targets -> functional_targets
    '''

    def test_pty(self):
        _ = self.signature_datasets
        _ = self.signature_datasets_meta
        _ = self.visualise_datasets
        _ = self.signature_heatmap
        self.export

    order = ['ZTIME','GTYPE','TEMP',]
    P = PTY = property
    # P = P
    from util import cached_property as CP
    import synotil.util as sutil 
    import pymisca.vis_util as pyvis
    plt = pyvis.plt    
    SRCDIR = os.path.realpath(os.path.dirname(__file__))

    

    ###############################
    #### Data loading #############
    inputs = Inputer() 
    # input_rnaseq  = None
    # input_datasets_meta = None
    # input_markers_df = None
    # input_rnaseq = "{input_rnaseq}"
    # input_datasets_meta = "{input_datasets_meta}"
    # input_markers_df = "{input_markers_df}"
    # FNAME = "{FNAME}"
    # OFNAME = "{OFNAME}"
    def __init__(self):
        _ = '''
        #### Parsing Command line arguments
        '''
        for x in "input_rnaseq, input_datasets_meta, input_markers_df, FNAME, OFNAME, TARGET".split(','):
            x = x.strip()
            k = "--%s"%x
            got = None
            if k in sys.argv: 
                got =  sys.argv[ sys.argv.index(k) + 1] 
            if getattr(self, x, got) is None:
                assert 0, "Missing --%s.\n Calling: %s"%(x,sys.argv)
            if got is not None:
                setattr(self, x, got)

        super(rnaseq_figure,self).__init__()

    # input_rnaseq = inputs._get_file('/home/feng/envs/Fig_POLYQ/rnaseq.pk')
    # input_datasets_meta = inputs._get_file('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv')
    # input_markers_df = inputs._get_file(pyext.f('{SRCDIR}/key_ath.csv'),raw=1)

    @CP
    def rnaseq(self):
        rnaseq = rnaseq_raw = pyext.readData(self.input_rnaseq)
        # rnaseq = rnaseq_raw = pyext.readData(pyext.f('{SRCDIR}/static.envs.Fig_POLYQ.rnaseq.pk'))
        rnaseq = rnaseq.copy()
        rnaseq.loc[:] = rnaseq.apply(pyext.log2p1)
    # rnaseq=rnaseq 
    # job['rnaseq'] = rnaseq
        return rnaseq

    @CP
    def datasets_meta(self):
        mcurr0 = pyext.readData(self.input_datasets_meta)
        # mcurr0 = pyext.readData( pyext.f('{SRCDIR}/static.results.0318-makeRNA-polyQ.mcurr0.csv') )
        mcurr0.columns = mcurr0.columns.str.upper()
        mcurr0['DISP_NAME']= pyext.df__format(mcurr0,'{TEMP}-{ZTIME}-{GTYPE}-{index}')
        return mcurr0
    # job['datasets_meta'] = mcurr0

    @CP
    def markers_df(self):
        keyDF = pyext.readData(self.input_markers_df)
        markers = ['LUX']
        markers_df = keyDFC = keyDF.query('BioName in %s'%markers)
        return markers_df
    # job['markers_df'] = markers_df = keyDFC = keyDF.query('BioName in %s'%markers)


    @property
    def signature_datasets(self, ):
        # datasets_meta, order):
    #     keys = ['ZTIME','GTYPE','TEMP',]
        keys = self.order
        mcurr0 = self.datasets_meta
        mcurr = mcurr0
        mcurr = mcurr.query("TEMP in ['22C','27C']")
        mcurr = mcurr.query("RUNID in ['181R','184R','177R']")
        mcurr = mcurr.query("ZTIME in ['ZT12']")
        mcurr = mcurr.query("GTYPE in ['COL','ELF3','ELF3-OX AT']")
        mcurr = mcurr.query("~BNAME.str.contains('adult')")
        mcurr = mcurr.sort_values(keys,)
        # mcurr.index.isin(rnaseq.columns)
        print(mcurr.shape)
        mcurr[keys]

        signature_datasets = self.rnaseq.reindex(columns = mcurr.index)
        return signature_datasets

    @property
    def signature_datasets_meta(self):
        mcurr = self.datasets_meta.reindex( self.signature_datasets.columns)
        return mcurr


    @PTY
    def visualise_datasets(self):
        self.to_override["visualise_datasets"] = None
        keys = self.order
        mcurr0 = self.datasets_meta
        mcurr = mcurr0
        mcurr = mcurr.query("RUNID in '181R,184R,177R,178R'.split(',')")
        # mcurr = mcurr.query("TEMP in ['22C','27C']")
        mcurr = mcurr.sort_values(keys)

        visualise_datasets = rnaseq.reindex(columns=mcurr.index)
        csv = mcurr0[keys].reindex(visualise_datasets.columns)
        csv["IN_HEATMAP"] = 1
        csv["IN_BOXPLOT"] = 1
        csv.to_csv( "visualise_datasets.csv")
        # print csv.shape
        return visualise_datasets
    @PTY
    def signature_heatmap(self,):
        # key ,signature_datasets, signature_targets):
        vdf = self.signature_datasets.reindex( self.signature_targets)
        # fig
        vdf.loc[:] = self.sutil.meanNorm(vdf)
        vdf = vdf.T
        im  = self.pyvis.heatmap(vdf,figsize=[12,7])
        return self.plt.gcf()
    @PTY
    def plot_signature_profile(self):

        mcurr = self.datasets_meta.reindex( self.signature_datasets.columns)
        ax = self.plt.gca()
        vdf = self.rnaseq.reindex( self.markers_df.index, columns=mcurr.index)
        vdf.columns = mcurr['DISP_NAME']
        vdf = vdf.T
        # vdf.co
        vdf.plot(xticks=range(len(vdf)),rot='vertical',ax=ax)
        ax.grid(1)
        ax.set_xlabel('dataset')
        ax.set_ylabel('log2(1+TPM_or_CPM)')
        return ax



    signature_CUTOFF = 0.99
    @P
    def signature_targets(self):
        plt = self.plt
        pyvis = self.pyvis
        vdf = self.signature_datasets
        signature_score = vdf.dot(self.signature_profile)
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
            pyvis.abline(x0=self.signature_CUTOFF)
    #     CUTOFF = 0.99 ### top 1% of 
        _targets = vdf.index[ppf > self.signature_CUTOFF]
    #     with pyext.getPathStack([WORKDIR,key],force=1) as stack:
        pyext.printlines(_targets,  pyext.f( '_temp-signature_targets-.it'))
        return _targets

    @P
    def signature_profile(self):
        vdf = self.signature_datasets
        vdf = vdf.reindex(index=self.markers_df.index)
        vdf.loc[:] = self.sutil.meanNorm(vdf)
        signature_profile = vdf.mean(axis=0)    
        return signature_profile


    @P
    def visualise_heatmap(self):
        vdf = self.visualise_datasets.reindex(
            self.signature_targets)
        
        vdf.columns = self.datasets_meta.reindex(vdf.columns)['DISP_NAME']
        vdf.loc[:] = self.sutil.meanNorm(vdf)
        vdf = vdf.T
        im = self.pyvis.df__heatmap(vdf,
                               cname='value',
                          )


        fig = self.plt.gcf()
        fig.set_figwidth(10)
        fig.set_figheight( 0.25*
            len(self.visualise_datasets.columns))
        fig.axes[0].set_title(self.FNAME)

        return fig

    @P
    def visualise_boxplot(self):

        vdf = self.visualise_datasets.reindex( 
            self.signature_targets)
        vdf.columns = self.datasets_meta.reindex(vdf.columns)['DISP_NAME']
        
        vdf.loc[:] = self.sutil.meanNorm( vdf )
    #     vdf = vdf.T
        fig,ax = self.plt.subplots(1,1,figsize=[16,7])
        vdf.boxplot(rot='vertical',ax=ax)
        ax.set_ylabel('meanNorm(log2( 1+ CPM_or_TPM))')

        fig.set_figwidth(7)
        fig.set_figheight(7)
        fig.axes[0].set_title(self.FNAME)

        return fig

    @P
    def export(self):
        df = pd.DataFrame({"rnaseq_signature_targets": self.signature_targets})
        with pyext.getPathStack([self.OUTDIR],force=1):
            df.to_csv("gene_lists_dataframe.csv",index=0)
        pass

    def main(self):
        import sys
        # if "--print-inputs" in sys.args:
        #     for x in self.inputs:
        #         print(x)
        #     return 0
        if "--run" in sys.argv:
            self.run()
            return 0
        else:
            assert 0,sys.argv




    @property   
    def visualise_datasets(self):
        df  = pyext.readData(self.FNAME).dropna()
        return self.rnaseq.reindex(columns=df.index)
    def run(self):
        fig = getattr(self, self.TARGET);


        # fig.set_figwidth(7)
        # fig.set_figheight(7)
        # fig.axes[0].set_title(self.FNAME)
        assert self.is_ready(), self.to_override

        fig__saveFile(fig, self.OFNAME)

if __name__ == '__main__':
    rnaseq_figure().main()
            # return json.dumps()
if __name__ =='__main__':
    if 0:

        job = rnaseq_figure()
        job.OUTDIR = "test-out"
        job.test_pty()
        job.plot_signature_profile
        job.signature_targets
        job.visualise_boxplot
        print(job.signature_datasets.head())

        import sys
        sys.exit(0)
