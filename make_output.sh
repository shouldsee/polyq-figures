rm -rf OUTPUT
mkdir -p OUTPUT
set -e
python2 src/make_boxplot.py
python2 src/make_heatmap.py
python2 src/make_chipseq_targets.py
python2 src/make_chipseq_pileups.py
python2 src/make_gene_lists.py
mkdir -p _build_temp
mv -f chip* home* genes* infiles* job* temp* _temp* _build_temp || true