from lazydict import LazyDictionary
import pymisca.ext as pyext

obj_chipseq_db = job = LazyDictionary()
job['TYPE'] = 'CHIPSEQ_DB'
job['SOURCE_DF'] = lambda s,k:[[][0],'must specify SOURCE_DF']
@pyext.setItem(job,'_SOURCE_DF')
def _func(job,key, SOURCE_DF):
    res = SOURCE_DF.copy()
    res.columns = res.columns.str.replace('[^a-zA-Z0-9_]','_').str.upper()
    assert "DATA_ACC" in res.columns,(res.columns,)
    return res

@pyext.setItem(job,"init")
def _func(job,key, _SOURCE_DF):
#     assert "DATA_ACC" in SOURCE_DF.columns,(SOURCE_DF.columns)

    @pyext.setItem(job)
    @pyext.PlainFunction
    def get__chipseq__bwfile( DATA_ACC):
        return get__chipseq__meta(DATA_ACC)['RPKMFILE']
    
    @pyext.setItem(job)
    @pyext.PlainFunction
    def get__chipseq__meta( DATA_ACC, SOURCE_DF = _SOURCE_DF):
        res = SOURCE_DF.set_index('DATA_ACC').loc[DATA_ACC].to_dict()
        return res
@pyext.setItem(job,'from/file')
# def _func(job,key,FNAME)
@pyext.PlainFunction
def from__file(FNAME,job=job):
    job = job.copy()
    df = pyext.readData('/home/feng/meta/meta_chip.tsv')
    df['DATA_ACC'] = df.index
    job['SOURCE_DF'] = df
    job['init']
    job['get__chipseq__bwfile']('192CS17')    
    return job

    