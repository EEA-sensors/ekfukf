#!/bin/sh

if test "$1" = ""; then
  echo "Usage: ./ekfukf_pack.sh x_y"
else

  echo Making temporary copy ...

  mkdir ekfukf
  cp -r src/* ekfukf

  echo Generating ekfukf_$1.tar.gz ...
  tar -cf ekfukf_$1.tar ekfukf/*
  gzip ekfukf_$1.tar

  echo Generating ekfukf_$1.zip ...
  zip -r ekfukf_$1.zip ekfukf/*
  rm -rf ekfukf

  echo Deleting temporary copy ...
  echo Done.
fi
