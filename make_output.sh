DIR=${1:-$PWD/OUTPUT}
rm -rf $DIR
rm -rf docs

mkdir -p $DIR
set -e
python2 src/make_boxplot.py  ## lucked
python2 src/make_heatmap.py  ## lucked
python2 src/make_chipseq_targets.py ##lucked
python2 src/make_chipseq_pileups.py ##lucked
python2 src/make_gene_lists.py ## lucked
python2 src/make_dep_graph.py
# mkdir -p _build_temp
mv -f chip* home* genes* infiles* job* temp* _temp* _build_temp || true

cp -prf src/ -t $DIR; rm -rf $DIR/src/static/
ln README.md $DIR

find $DIR -type l -delete  ### jekyll does not work with symlink
cd $DIR
URL=https://gist.githubusercontent.com/glowinthedark/b1f5900be2490c5371f827a49fd09f49/raw/db45a596db2274b94206f4dfa6d479cee4e49845/generate_directory_index.py
curl -L $URL | python2 -
cd $OLDPWD

mv $DIR docs
echo "[DONE]"
