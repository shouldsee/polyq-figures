# %%python2
# from rnaseq_figure import job as template
from luck.shorts import RNS, TSSR, NCR, LSC
RULE = TSSR
ns = RNS()

from util import _get_output_file


# @cached_property
# !mkdir -p results
for FNAME in [
    'src/BOXPLOT_1_OX-lines.tsv',
    'src/BOXPLOT2_mutants.tsv',
    'src/BOXPLOT3_Qlines.tsv',
]:
    import os
    # job = template.copy()
    FNAME = FNAME 
    OFNAME = "OUTPUT/" + FNAME.replace('/','_')+'.svg'
    OFNAME = _get_output_file(OFNAME)
    TEMP_SCRIPT = OFNAME+'.py'
    import shutil
    # shutil.copy2("rnaseq_figure_2.py",TEMP_SCRIPT)
    with open("rnaseq_figure_2.py","r") as fold:
        with open(TEMP_SCRIPT,"w") as f:
            f.write(
f'''
import sys,os
os.chdir("{os.path.realpath(os.getcwd())}"); sys.path.append("");
{fold.read()}
# from rnaseq_figure_2 import rnaseq_figure as template
import pymisca.ext as pyext
# from util import _get_output_file


def fig__saveFile(fig, ofname, transparent = False, **kwargs):
    fig.savefig(ofname,
                bbox_inches='tight',
                        transparent= transparent,
                        facecolor=fig.get_facecolor(),
                        **kwargs)
    return 1

class rnaseq_figure(rnaseq_figure):
    @property   
    def visualise_datasets(self):
        df  = pyext.readData(self.FNAME).dropna()
        return self.rnaseq.reindex(columns=df.index)
    def run(self):
        # fig = job['visualise_boxplot'];
        fig = self.visualise_boxplot;
        fig.set_figwidth(7)
        fig.set_figheight(7)
        fig.axes[0].set_title(self.FNAME)
        assert self.is_ready(), self.to_override

        fig__saveFile(fig, self.OFNAME)

if __name__ == '__main__':
    rnaseq_figure.main()
    ''',
)

    # LSC("python2 {TEMP_SCRIPT} > {TEMP_SCRIPT}.temp")
    RULE.MWF(ns, OFNAME, [FNAME] + LSC(f"python2 {TEMP_SCRIPT} --print-inputs").decode().replace("\n"," "),  f"python2 {TEMP_SCRIPT} --run" )
    # list( self.inputs), self.run)
    break


#     fig;
# job['visualise_datasets'].head()