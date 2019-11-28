
from pymisca.atto_jobs import ModuleJob
import path
import os
from util import _get_file,_get_output_file
import src.util as _util

FORCE = 1
SRC_DIR = os.path.dirname(__file__)# 
    
# _SRC = path.Path(__file__).dirname()
RESULTS = res = [
    ModuleJob({
    "MODULE_FILE": "./src/plot_bigwig_pileup.py",
    "DATA":{
        "PEAK_FILE": _util._get_middle_file('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),

#         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
        "OUTDIR": _util._get_output_file("OUTPUT/fig-2c"),
        'FORCE':FORCE,
        "BIGWIG_FILES":[
            _get_file('Mapped_data/192C/S17/ELF3myc-17C-ZT10_S17_Ath-TAIR10_RPKM.bw'),
            _get_file('Mapped_data/192C/S18/ELF3myc-27C-ZT10_S18_Ath-TAIR10_RPKM.bw'),
            
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S17/ELF3myc-17C-ZT10_S17_Ath-TAIR10_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S18/ELF3myc-27C-ZT10_S18_Ath-TAIR10_RPKM.bw'
        ],
        "AXIS_DICTS":[
            {'set_ylim':[0,12],
            'set_ylabel':"log2( RPKM_at_peak )"},
            {'set_xlabel':"distance to peak(bp)",
             'set_ylabel':"RPKM",
             'set_ylim':[0,800],
            }
        ],
}
}),
    ModuleJob({
    "MODULE_FILE":"./src/plot_bigwig_pileup.py",
    "DATA":{
        "PEAK_FILE": _util._get_middle_file('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),
        
#         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
        "OUTDIR": _util._get_output_file("OUTPUT/fig-2d"),
        'FORCE':FORCE,
        "BIGWIG_FILES":[
            _get_file('Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw'),
            _get_file('Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw'),
            _get_file('Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw'),
            _get_file('Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw'),
            
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw',
        ],

        "AXIS_DICTS":[
            {'set_ylim':[0,12],
            'set_ylabel':"log2( RPKM_at_peak )"},
            {'set_xlabel':"distance to peak(bp)",
             'set_ylabel':"RPKM",
             'set_ylim':[0,800],
            }
        ],
}
}),
    ModuleJob({
    "MODULE_FILE":"./src/plot_bigwig_pileup.py",
    "DATA":{
        "PEAK_FILE": _util._get_middle_file('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),
#         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
        "OUTDIR": _util._get_output_file("OUTPUT/fig-2e"),
        'FORCE':FORCE,
        "BIGWIG_FILES":[
            _get_file('Mapped_data/189C/S10/gELF3myc-17C_S10_RPKM.bw'),
            _get_file('Mapped_data/189C/S11/gELF3myc-27C_S11_RPKM.bw'),
            _get_file('Mapped_data/189C/S16/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-17C_S16_RPKM.bw'),
            _get_file('Mapped_data/189C/S17/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-27C_S17_RPKM.bw'),                        
            
#             _getFile("189C/S10/gELF3myc-17C_S10_RPKM.bw"),
           
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S10/gELF3myc-17C_S10_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S11/gELF3myc-27C_S11_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S16/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-17C_S16_RPKM.bw',
#             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S17/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-27C_S17_RPKM.bw'
        ],
            
        "AXIS_DICTS":[
            {'set_ylim':[0,12],
            'set_ylabel':"log2( RPKM_at_peak )"},
            {'set_xlabel':"distance to peak(bp)",
             'set_ylabel':"RPKM",
             'set_ylim':[0,800],
            }
        ],
}
}),
    
]