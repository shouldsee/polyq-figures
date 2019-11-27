from rnaseq_figure import job as template
import pymisca.ext as pyext
import src.util as _util

def fig__saveFile(fig, ofname, transparent = False, **kwargs):
    fig.savefig(ofname,
                bbox_inches='tight',
                        transparent= transparent,
                        facecolor=fig.get_facecolor(),
                        **kwargs)
    return 1
for FNAME in [
    'src/FIG-S6_HEATMAP1_Q-lines.tsv',
    'src/FIG-S5_HEATMAP2_genotypes.tsv',
    'src/FIG-S7_HEATMAP3.tsv',
#     'src/BOXPLOT_1_OX-lines.tsv',
#     'src/BOXPLOT2_mutants.tsv',
#     'src/BOXPLOT3_Qlines.tsv',
    
]:
    print( FNAME)
    job = template.copy()
    job['WORKDIR'] = "./WORKDIR"

    @pyext.setItem(job,'visualise_datasets')
    def _func(self,key,datasets_meta,rnaseq):
        df  = pyext.readData(FNAME).dropna()
        return rnaseq.reindex(columns=df.index[::-1])

    fig = job['visualise_heatmap'];
    fig.set_figwidth(10)
    fig.set_figheight( 0.25*
        len(job['visualise_datasets'].columns))
    fig.axes[0].set_title(FNAME)
    
    with pyext.getPathStack(['OUTPUT'],force=1):
        OFNAME = FNAME.replace('/','_')+'.png'
        OFNAME = _util._get_output_file(OFNAME)
        fig__saveFile(fig,OFNAME)
    fig;
    
#     fig.savefig('results/'+FNAME.replace('/','_')+'.png')
#     fig;