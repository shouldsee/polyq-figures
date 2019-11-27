'''
Produce Soft files for editing

'''
import pymisca.ext as pyext
node = pyext.relFile__asModule('./get_meta_soft.py')
dfc = node.df_figureRegistry()

# import pymisca.ext as pyext
# node = pyext.file__asModule('src/get_meta_soft.py')
node.WORKDIR = lambda: pyext.path.Path('/home/feng/envs/0726-polyq/WORKDIR.submit/').realpath()
node.get_soft_text()
pyext.MDFile('./WORKDIR.submit/get_soft_text.tar.gz')