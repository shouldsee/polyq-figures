import path
import os

# import pymisca.header
# from filelock import FileLock
import json
from collections import OrderedDict as _dict
import collections
import time
SRC_DIR = os.path.dirname(__file__)

STATIC_DIR = 'src/static'
INDEX_FILE = os.path.realpath("OUTPUT/index.json.list")


class Inputer(object):
    def __init__(self, ):
        self.input_files = {}
        return 

    def __repr__(self):
        return repr(self.input_files)
    def __iter__(self):
        return iter(self.input_files.keys())

    def _get_file(self,fn, *a,**k):
        fn = _get_file(fn,*a,**k)
        self.input_files[fn]= None
        return fn

# import pymisca.ext as pyext
# from pymisca.events import LinkEvent,CopyEvent
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


def _get_output_file(FNAME,INDEX_FILE=INDEX_FILE):
    d = [
        ("TS",time.time()),
        ("OUTPUT_FILE",FNAME,),
        # ("RUNTIME_FILE", pymisca.header.name__lookup('__file__',level=-1)),
    ]
    d = _dict(d)
    # print (json.dumps(d,indent=4))
#     with FileLock(  INDEX_FILE +'.lock'):
#         with open(INDEX_FILE, "a") as f:
# #             f.write(json.dumps(d)+'\n')
#             json.dump(d, f)
#             f.write('\n')
    return FNAME
get_output_file = _get_output_file

def _get_file(FNAME, raw =False):
    if ON_SERVER or raw:
#         FNAME = '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/%s'%
        FNAME = FNAME
    else:
        # FNAME = path.Path(SRC_DIR)/ "static" / FNAME.lstrip('/').replace('/','.') ##### caution_1
        FNAME = path.Path(SRC_DIR)/ "static" / FNAME.lstrip('/').replace('/','.') ##### caution_1
        # FNAME = FNAME.relpath(os.getcwd())
        # print('[dbg]',SRC_DIR, FNAME)
    assert os.path.exists(FNAME),FNAME
    return str(FNAME)

def _get_middle_file(FNAME):
#     FNAME = 'OUTPUT/'+FNAME
    assert os.path.exists(FNAME),FNAME
    return FNAME



class cached_property(object):
    """
    Descriptor (non-data) for building an attribute on-demand on first use.
    """
    def __init__(self, factory):
        """
        <factory> is called such: factory(instance) to build the attribute.
        """
        self._attr_name = factory.__name__
        self._factory = factory

    def __get__(self, instance, owner):
        # Build the attribute.
        attr = self._factory(instance)

        # Cache the value; hide ourselves.
        setattr(instance, self._attr_name, attr)

        return attr

