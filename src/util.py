import path
import os

SRC_DIR = os.path.dirname(__file__)

STATIC_DIR = 'src/static'

import os
# import pymisca.ext as pyext
from pymisca.events import LinkEvent,CopyEvent
def _file_static_copy( IN_DIR, IN_BASENAME, OUT_DIR, OUT_BASENAME):
    if IN_BASENAME is None:
        IN_BASENAME = IN_DIR
    if OUT_BASENAME is None:
        OUT_BASENAME = '.'.join(IN_BASENAME.lstrip('/').split('/'))  #### caution_1
    return CopyEvent(
        os.path.join(IN_DIR,IN_BASENAME),
        os.path.join(OUT_DIR,OUT_BASENAME),
        force=1)

# ON_SERVER = 1
ON_SERVER = 0
def _get_file(FNAME):
    if ON_SERVER:
#         FNAME = '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/%s'%
        FNAME = FNAME
    else:
        FNAME = path.Path(SRC_DIR)/ "static" / FNAME.lstrip('/').replace('/','.') ##### caution_1
    assert os.path.exists(FNAME),FNAME
    ###
    return FNAME

def _get_output_file(FNAME):
    FNAME = 'OUTPUT/'+FNAME
    assert os.path.exists(FNAME),FNAME
    return FNAME
