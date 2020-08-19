from pymisca.plotters import plotters
# from src.obj_chipeq_db import obj_chipseq_db
import synotil.jobs as sjob
from pymisca.atto_jobs import AttoJob

import pymisca.ext as pyext

class plot_bigwig_pileup(AttoJob):
    PARAMS_TRACED = [
        ('BIGWIG_FILES',('list:AttoPath',[])),
        ('BIGWIG_NAMES',('list:unicode',[])),
        ('BIGWIG_PAIRS',('list:list:object',[])),
        
        ('PEAK_FILE',('AttoPath','')),
        ('OUTDIR',('AttoPath','')),
        
        ('NCORE',('int',4)),
        ('FORCE',('int',0)),
        
        ('AXIS_DICTS',('list:dict:object',[])),
        ('PARAMS',('dict:object',{
            'outerRadius':500,
        })),
    ]
    def _run(self):
        kw = self._data
        BIGWIG_PAIRS = kw['BIGWIG_PAIRS']
        BIGWIG_FILES = kw['BIGWIG_FILES']
        BIGWIG_NAMES = kw['BIGWIG_NAMES']
        PARAMS = kw['PARAMS']
        FORCE = kw['FORCE']
        NCORE = kw['NCORE']
#         PARAMS['NCORE'] = kw['NCORE']
#         kw['LAST_FILE'] = kw['OFNAME'] = OFNAME = kw['OFNAME'].realpath()
#         kw['LAST_DIR'] = OFNAME.dirname()
        kw['PEAK_FILE'] = PEAK_FILE = kw['PEAK_FILE'].realpath()
        assert kw['OUTDIR']
        kw['OUTDIR'] = OUTDIR = kw['OUTDIR'].realpath()
        AXIS_DICTS = kw['AXIS_DICTS'] 
        
        if BIGWIG_PAIRS:
            assert not BIGWIG_NAMES and not BIGWIG_FILES,('conflict arguments')
            BIGWIG_FILES, BIGWIG_NAMES = zip(BIGWIG_PAIRS)
            pass
        else:
            assert BIGWIG_FILES
            if not BIGWIG_NAMES:
                BIGWIG_NAMES = map(pyext.getBname,BIGWIG_FILES)
                
        with pyext.getPathStack([OUTDIR],force=1):
            OFNAME = 'DONE'
    #         OFNAME = 
            if not FORCE and pyext.file__notEmpty(OFNAME):
                pass
            else:
                figs,(bwTable,bwAvg) = sjob.figs__peakBW(
                    peakFile=PEAK_FILE,
                    bwFiles = BIGWIG_FILES,
                    outIndex = None if not len(BIGWIG_NAMES) else BIGWIG_NAMES,
                    NCORE=NCORE,**PARAMS )
                fig = figs.values()[0]
                for (axis,d_ax) in zip( fig.axes, AXIS_DICTS):
                    pyext.obj__dict__call(axis,d_ax)

                    
                for OFNAME in ['MAIN.png','MAIN.svg']:
                    plotters.fig__save(fig,OFNAME)
                bwTable.to_pickle('bwTable.pk')
                bwAvg.to_pickle('bwAvg.pk')
                pyext.printlines([pyext.dppJson(kw)],"DATA.json")
                pyext.printlines(['DONE'],'DONE')

#         if not OFNAME:
#         if not OFNAME
#     "outerRadius":500,
#     "NCORE":4,
#     "peakFile":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
#     "bwFiles":lambda self,key,DATA_ACC_LIST,chipseq_db:
#         map(chipseq_db['get__chipseq__bwfile'],DATA_ACC_LIST)
#     "OFNAME":
        
#     return 