# polyq-figures

Results: https://shouldsee.github.io/polyq-figures/

## Running the code

- clone the repo to local `git clone https://github.com/shouldsee/polyq-figures && cd polyq-figures`
- Make sure the dependencies were installed through
```
sudo apt-get install -y graphviz
sudo bash -e install-binary.sh
pip install requirements.txt --user
```




- [New] Install `wget https://github.com/shouldsee/luckmake/releases/download/0.0.6/luckmake && chmod +x luckmake`luckmake and run `./luckmake OUTPUT`
- [Old]: run the main script `bash make_output.sh`
- inspecting the outputed folder `./OUTPUT/`

## Overview

This repository contains code to reproduce figures in "A prion-like domain in ELF3 functions as a thermosensor in Arabidopsis" (in press)

You can view the results at https://shouldsee.github.io/polyq-figures/ 


```
src/static  # contains processed data
src/        # contains scripts to produce figues
docs/       # contains output for static sites
```

# CHANGELOG

- 20200715: reindexed pileup figures. added "fig2e-inlet"

# Dependency Between tasks

`luckmake graph`

![graph](LUCKFILE.py.dot.svg)
