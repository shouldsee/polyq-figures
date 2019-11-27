shopt -s extglob ###expanding !($HASH)
### password is reviewpaw
DIR=${1:-$PWD}
HASH=a93194bed5709a448522d69b9a722c02d1e7e71f
find $DIR -type l -delete  ### jekyll does not work with symlink
mkdir -p $HASH && mv ./!($HASH) ./$HASH
curl -L https://raw.githubusercontent.com/matteobrusa/Password-protection-for-static-pages/master/index.html >index.html  

echo SITE encrypted
