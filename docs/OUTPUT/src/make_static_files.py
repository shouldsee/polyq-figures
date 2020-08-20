
from src.util import _file_static_copy,STATIC_DIR    
it = [
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/192C/S17/ELF3myc-17C-ZT10_S17_Ath-TAIR10_RPKM.bw',
     STATIC_DIR,
    None),
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/192C/S18/ELF3myc-27C-ZT10_S18_Ath-TAIR10_RPKM.bw',
     STATIC_DIR,
    None),
    
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw',
     STATIC_DIR,
    None),   
    
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw',
     STATIC_DIR,
    None),   
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw',
     STATIC_DIR,
    None),   
    
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw',
     STATIC_DIR,
    None),   
    

    
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/189C/S10/gELF3myc-17C_S10_RPKM.bw',
     STATIC_DIR,
    None),   
    
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/189C/S11/gELF3myc-27C_S11_RPKM.bw',
     STATIC_DIR,
    None),   
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/189C/S16/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-17C_S16_RPKM.bw',
     STATIC_DIR,
    None),   
    ('/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/',
     'Mapped_data/189C/S17/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-27C_S17_RPKM.bw',
     STATIC_DIR,
    None),       
    
]

[ _file_static_copy(*x) for x in it]
# print 
"DONE"



######
import src.get_meta_soft as ns
DATA_ACC_LIST = [
    "189CS10",
    "189CS11",
    "192CS17",
    "192CS18"]

meta = ns.df_mappedData_chipseq()
# bwFile = meta.loc[DATA_ACC,'bw']
bwFiles = meta['bw']
bwFiles=bwFiles[ DATA_ACC_LIST].tolist()
npkFile = meta.loc["192CS17",'narrowPeak']

it = [(x,None,STATIC_DIR,None) for x in bwFiles+[npkFile]]
[ _file_static_copy(*x) for x in it]
"DONE"
# map(lambda x:file_static_link(*x),it)

'''
cat src/get_meta_soft.py | grep -n _readData
21:_readData = lambda x,**kw: pyext.readData(_get_file(x),**kw)
24:#         _readData("/home/feng/envs/0726-polyq/src/BOXPLOT_1_OX-lines.tsv").dropna().index.tolist(),
25:#         _readData("/home/feng/envs/0726-polyq/src/BOXPLOT2_mutants.tsv").dropna().index.tolist(),
26:#         _readData("/home/feng/envs/0726-polyq/src/BOXPLOT3_Qlines.tsv").dropna().index.tolist(),
39:    mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
46:    mcurr0 = _readData('/home/feng/meta/meta_chip.tsv',guess_index=0)
185:    mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
224:    mcurr0 = _readData('/home/feng/meta/meta_chip.tsv',guess_index=0)
225:#     mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
348:    mcurr = _readData(
353:#     mcurr = _readData(
367:    df = _readData('/home/feng/envs/upGeo/results/0407-database-mapped/mcurr.csv')
653:        df = _readData(sample['npkfile'],columns=pyext.columns.bed,header=None)
666:#     mcurr = _readData(
678:    mappedMeta = _readData(
'''

it = [
  #### metafile databases specific to slcu server
 '/home/feng/envs/upGeo/results/0407-database-mapped/mcurr.csv',
 '/home/feng/envs/upGeo/results/0424-database-raw/mcurr.csv',
 '/home/feng/meta/meta_chip.tsv',
 '/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',
    
  ### Genome reference
  "/home/feng/ref/Arabidopsis_thaliana_TAIR10/Annotation/genes.gtf.cds",
  "/home/feng/ref/Arabidopsis_thaliana_TAIR10/genome.sizes",
    
  #### rnaseq files
  '/home/feng/envs/Fig_POLYQ/rnaseq.pk',
    
]

from src.util import _file_static_copy,STATIC_DIR
it = [(x,None,STATIC_DIR,None) for x in it]
[ _file_static_copy(*x) for x in it]