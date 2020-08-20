
import pymisca.ext as pyext
pyext.os.chdir('/home/feng/envs/0726-polyq/')
ns = pyext.file__asModule('/home/feng/envs/0726-polyq/src/get_meta_soft.py')
WORKDIR = ns.WORKDIR
OUTDIR = WORKDIR() / "get_soft_text"

pyext.real__dir(dirname = OUTDIR)
_samples = ns.sample_init_full()
for sample in _samples:
    print '[template_finalise]',sample['data_acc'],'...'
    ns.sample_template_finalise(sample)
    res = ns.sample_get_depth(sample, WORKDIR()/sample["data_acc"]/"sample_get_depth.json")
    sample.update(res )
    break
    
