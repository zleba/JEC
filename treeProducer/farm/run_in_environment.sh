#!/bin/zsh

export LD_LIBRARY_PATH=$NAFTMP:$LD_LIBRARY_PATH_STORED
export PATH=$NAFTMP:$PATH

cd $NAFTMP

echo $*

$*

echo Done
