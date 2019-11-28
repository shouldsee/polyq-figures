
from distutils.core import setup

setup(name='polyq_figures',
      version='0.0.1',
      description='Scripts for reproducing polyQ analysis',
      author='Feng Geng',
      author_email='fg368@cam.ac.uk',
      url='none',
      install_requires=[ x.strip() for x in open("requirements.txt","r") 
        if x.strip() and not x.strip().startswith("#") ],
      package_dir = {},
      packages=[],
     )
