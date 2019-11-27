DIR_BIN=/usr/local/bin
URL=https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
wget -nc $URL
chmod +x bedtools; ln -sft $DIR_BIN $PWD/bedtools
