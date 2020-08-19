# %%python2
# from rnaseq_figure import job as template
from luck.shorts import RNS, TSSR, NCR, LSC
RULE = TSSR
ns = RNS()

from util import _get_output_file, Inputer
import os
SRCDIR = os.path.realpath(os.getcwd())


if 1:
    inputs = Inputer()
    input_rnaseq = inputs._get_file('/home/feng/envs/Fig_POLYQ/rnaseq.pk')
    input_datasets_meta = inputs._get_file('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv')
    input_markers_df = inputs._get_file(f'{SRCDIR}/src/key_ath.csv',raw=1)
    for FNAME in inputs:
        RULE.MWF(ns, FNAME, None)

    #### Making boxplots
    OUTS = []
    for FNAME in [
        'src/BOXPLOT_1_OX-lines.tsv',
        'src/BOXPLOT2_mutants.tsv',
        'src/BOXPLOT3_Qlines.tsv',
    ]:
        # job = template.copy()
        FNAME = FNAME 
        OFNAME = "OUTPUT/" + FNAME.replace('/','_')+'.svg'
        OFNAME = _get_output_file(OFNAME)
        SCRIPT = "src/rnaseq_figure_2.py"
        RULE.MWF(ns, FNAME, None)
        RULE.MWF(ns, OFNAME, ' '.join([FNAME] + [SCRIPT, input_rnaseq, input_datasets_meta, input_markers_df]),  f"\
            python2 {SCRIPT} \
            --input_rnaseq {input_rnaseq} \
            --input_datasets_meta {input_datasets_meta} \
            --input_markers_df {input_markers_df} \
            --FNAME {FNAME} \
            --OFNAME {OFNAME} \
            --TARGET visualise_boxplot \
            --run \
            ")
        OUTS.append(OFNAME)
    NCR.M(ns, "boxplots", ' '.join(OUTS))

    OUTS = []
    for FNAME in [
        'src/FIG-S6_HEATMAP1_Q-lines.tsv',
        'src/FIG-S5_HEATMAP2_genotypes.tsv',
        'src/FIG-S7_HEATMAP3.tsv',
    ]:    
        OFNAME = "OUTPUT/"+FNAME.replace('/','_')+'.png'
        OFNAME = _util._get_output_file(OFNAME)

        RULE.MWF(ns, OFNAME, ' '.join([FNAME] + [SCRIPT, input_rnaseq, input_datasets_meta, input_markers_df]),  f"\
            python2 {SCRIPT} \
            --input_rnaseq {input_rnaseq} \
            --input_datasets_meta {input_datasets_meta} \
            --input_markers_df {input_markers_df} \
            --FNAME {FNAME} \
            --OFNAME {OFNAME} \
            -- TARGET visualise_heatmap \
            --run \
            ")
        OUTS.append(OFNAME)
    NCR.M(ns, "heatmaps", ' '.join(OUTS))


NCR.M(ns, "all", "boxplots heatmaps", )

