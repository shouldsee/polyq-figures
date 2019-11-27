# %%python2
from rnaseq_figure import job as template
import pymisca.ext as pyext
def fig__saveFile(fig, ofname, transparent = False, **kwargs):
    fig.savefig(ofname,
                bbox_inches='tight',
                        transparent= transparent,
                        facecolor=fig.get_facecolor(),
                        **kwargs)
    return 1
# !mkdir -p results
for FNAME in [
    'src/BOXPLOT_1_OX-lines.tsv',
    'src/BOXPLOT2_mutants.tsv',
    'src/BOXPLOT3_Qlines.tsv',
]:
    job = template.copy()

    @pyext.setItem(job,'visualise_datasets')
    def _func(self,key,datasets_meta,rnaseq):
        df  = pyext.readData(FNAME).dropna()
        return rnaseq.reindex(columns=df.index)
    #     return df

    fig = job['visualise_boxplot'];
    fig.set_figwidth(7)
    fig.set_figheight(7)
    fig.axes[0].set_title(FNAME)
#     OFNAME  = 'OUTPUT/'+FNAME.replace('/','_')+'.svg'
    with pyext.getPathStack(['OUTPUT'],force=1):
        OFNAME = FNAME.replace('/','_')+'.svg'
        fig__saveFile(fig,OFNAME)
    fig;
# job['visualise_datasets'].head()