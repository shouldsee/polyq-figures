shopt -s extglob ###expanding ./!($HASH)
### password is reviewPAW
DIR=${1:-$PWD}
HASH=a93194bed5709a448522d69b9a722c02d1e7e71f
find $DIR -type l -delete  ### jekyll does not work with symlink

### add index
URL=https://gist.githubusercontent.com/glowinthedark/b1f5900be2490c5371f827a49fd09f49/raw/db45a596db2274b94206f4dfa6d479cee4e49845/generate_directory_index.py
curl -L $URL | python -
mkdir -p $HASH && mv ./!($HASH) ./$HASH
curl -L https://raw.githubusercontent.com/matteobrusa/Password-protection-for-static-pages/master/index.html >index.html  


echo SITE encrypted
