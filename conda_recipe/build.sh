#!/bin/sh

mkdir -p $PREFIX/bin/sub_bin

make CXX=${CXX}
cp metaplatanus $PREFIX/bin
cp src/scripts/tgsgapcloser_mod $PREFIX/bin
cp -r sub_bin/* $PREFIX/bin/sub_bin

wget https://github.com/rkajitani/NextPolish/releases/download/v1.3.1/NextPolish.tgz
tar xf NextPolish.tgz
make -C NextPolish
cp -r NextPolish/nextPolish NextPolish/lib NextPolish/bin $PREFIX/bin/sub_bin
