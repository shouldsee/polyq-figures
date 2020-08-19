# %matplotlib inline
import pymisca.ext as pyext
import matplotlib
matplotlib.use('Agg')
import os
import shutil



from pymisca.util import colGroupMean
pyext.colGroupMean = colGroupMean
import pymisca.vis_util as pyvis

import src.util as _util
from util import _get_file
import synotil.dio
import synotil.qcplots

import sys

SRC_DIR = os.path.dirname(__file__)
# pyext.os.chdir('/home/feng/envs/0726-polyq/')
# ns = pyext.file__asModule('/home/feng/envs/0726-polyq/src/get_meta_soft.py')
ns = pyext.file__asModule( SRC_DIR + '/get_meta_soft.py')
meta = ns.df_mappedData_chipseq()
# DATA_ACC = "189CS11"
# npkFile = meta.loc["189CS10",'narrowPeak']
# npkFile = meta.loc["189CS10",'narrowPeak']

inputs = []
def _get_file(fn):
  fn= _util._get_file(fn)
  inputs.append(fn)
  return fn
outputs  = []
def _get_output_file(fn):
  fn = _util._get_output_file(fn)
  outputs.append(fn)
  return fn

_ = '''Get input files'''
npkFile = meta.loc["192CS17",'narrowPeak']
npkFile = _get_file(npkFile)
GSIZE = _get_file("/home/feng/ref/Arabidopsis_thaliana_TAIR10/genome.sizes")
cds_file = _get_file("/home/feng/ref/Arabidopsis_thaliana_TAIR10/Annotation/genes.gtf.cds")
DATA_ACC_LIST = [
    "189CS10",
    "189CS11",
    "192CS17",
    "192CS18"]
meta = ns.df_mappedData_chipseq()
bwFiles = meta['bw']
bwFiles = bwFiles[ DATA_ACC_LIST]
bwFiles = map(_get_file, bwFiles)
if '--print-inputs' in sys.argv:
  for x in inputs: print(x)
  sys.exit(0)

_ = '''Get output files'''
OUTPUT_BED_FILE = _get_output_file('OUTPUT/chipseq_differential_binding.peak_list.bed')
outputs.append(OUTPUT_BED_FILE+'.summit')
OUTPUT_CSV_FILE = _get_output_file("OUTPUT/chipseq_differential_binding.peak_list.csv")
OUTPUT_CSV_GENE_FILE = _get_output_file("OUTPUT/chipseq_targets_genes_job.peak_list.csv")

if '--print-outputs' in sys.argv:
  for x in outputs: print(x)
  sys.exit(0)

if '--run' not in sys.argv:
  sys.exit(1)


df = pyext.readData(npkFile,'tsv',header=None,columns=pyext.columns.bed)
df['FC'].apply(pyext.np.log2).hist(bins=30)
sel = df['neglogPval'] > 3.
print(pyext.np.sum(sel))
df = df.loc[sel].to_csv('temp.bed',sep='\t',header=None)



[bwFiles]
res = synotil.dio.extract_bigwig_multiple(bedFile='temp.bed',
                                    outIndex=DATA_ACC_LIST,
                                    bwFiles=bwFiles,
                                    radius=300,stepSize=10)




tab = colGroupMean(res)
tab = tab.apply(pyext.log2p1)
xs =(tab['189CS10']) - (tab['189CS11'])
ys = tab['192CS17'] - tab['192CS18']
clu = (xs + ys) > 2
pyvis.qc_2var( xs,ys,nMax=-1,xlim=[-2,4],ylim=[-2,4],clu=clu)

df = pyext.readData('temp.bed',columns=pyext.columns.bed,guess_index=0).set_index('acc',drop=0)

# OUTPUT_FILE = 'OUTPUT/0918-elf3target.bed'
pyext.dir__real(dirname="OUTPUT")

[OUTPUT_BED_FILE]
df.loc[clu].to_csv(OUTPUT_BED_FILE,sep='\t',header=None,index=0)

[OUTPUT_CSV_FILE]
df.loc[clu].to_csv( OUTPUT_CSV_FILE
                   ,sep=',',header=None,index=0)


#### making gene list
FNAME =  OUTPUT_BED_FILE
shutil.move(synotil.dio.npk_expandSummit(FNAME,radius=1,), FNAME+'.summit')
INPUT_BED_FILE = FNAME +'.summit'

# INPUT_BED_FILE = 'OUTPUT/chipseq_differential_binding.peak_list.bed.summit'


[GSIZE, cds_file]
cds_summit = synotil.dio.bed__leftSummit(cds_file,GSIZE=GSIZE)

res = synotil.qcplots.qc_summitDist(peak1= INPUT_BED_FILE,
                              CUTOFF=500,
                              peak2 = cds_summit,
                              GSIZE =GSIZE,
                             );

res[0].to_csv(OUTPUT_CSV_GENE_FILE)


# # OUTPUT_BED_FILE = os.path.realpath(OUTPUT_BED_FILE)
# with pyext.TempDirScope(keep=0) as d:
#     _TEMP_OUTPUT=  synotil.dio.job__nearAUG(peakFile = d._root / OUTPUT_BED_FILE,
#                   featFile = _get_file("/home/feng/ref/Arabidopsis_thaliana_TAIR10/Annotation/genes.gtf.cds"),
# #                  featFile = _get_file("/home/ref_genew/Arabidopsis_thaliana_TAIR10/Annotation/Genes/genes.gtf"),
# #                   featFile = _get_file("/home/ref_genew/Arabidopsis_thaliana_TAIR10/Annotation/Genes/genes.gtf"),
#                   GSIZE = _get_file("/home/feng/ref/Arabidopsis_thaliana_TAIR10/genome.sizes"),
#                   CUTOFF=500,
#                  )
#     ### csv to tsv
#     pyext.readData(_TEMP_OUTPUT).to_csv(d._root / "OUTPUT/chipseq_targets_genes_job.peak_list.csv",)
# #     shutil.move(_TEMP_OUTPUT, d._root / "OUTPUT/chipseq_targets_genes_job.peak_list.tsv")
# #     pyext.readData("OUTPUT/chipseq_targets_genes_job.peak_list.tsv")

# # OUTPUT_FILE = 
