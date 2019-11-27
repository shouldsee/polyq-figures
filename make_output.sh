rm -rf OUTPUT
python2 src/make_boxplot.py
python2 src/make_heatmap.py
python2 src/make_chipseq_targets.py
python2 src/make_chipseq_pileups.py
python2 src/make_gene_lists.py
rm -f chip* home* genes* infiles* job* temp* _temp*