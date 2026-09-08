#!/bin/bash

# Check if both hts_path and deflate_path are provided
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <hts_path> <deflate_path>"
  exit 1
fi

HTS_PATH=$1
DEFLATE_PATH=$2

# Update Makefile
sed -i "s|HTSLIB_INSTALL_PATH = .*|HTSLIB_INSTALL_PATH = $HTS_PATH|g" Makefile
sed -i "s|LIBDEFLATE_INSTALL_PATH = .*|LIBDEFLATE_INSTALL_PATH = $DEFLATE_PATH|g" Makefile

# Update env.sh
sed -i "s|HTSLIB_INSTALL_PATH=.*|HTSLIB_INSTALL_PATH=$HTS_PATH|g" env.sh
sed -i "s|LIBDEFLATE_INSTALL_PATH=.*|LIBDEFLATE_INSTALL_PATH=$DEFLATE_PATH|g" env.sh

echo "Updated HTSLIB_INSTALL_PATH in Makefile to $HTS_PATH"
echo "Updated LIBDEFLATE_INSTALL_PATH in Makefile to $DEFLATE_PATH"
echo "Updated HTSLIB_INSTALL_PATH in env.sh to $HTS_PATH"
echo "Updated LIBDEFLATE_INSTALL_PATH in env.sh to $DEFLATE_PATH"

