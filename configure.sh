#!/bin/bash

# Check if hts_path is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <hts_path>"
  exit 1
fi

HTS_PATH=$1

# Update Makefile
sed -i "s|HTSLIB_INSTALL_PATH = .*|HTSLIB_INSTALL_PATH = $HTS_PATH|g" Makefile

# Update env.sh
sed -i "s|HTSLIB_PATH=.*|HTSLIB_PATH=$HTS_PATH|g" env.sh

echo "Updated HTSLIB_INSTALL_PATH in Makefile to $HTS_PATH"
echo "Updated HTSLIB_PATH in env.sh to $HTS_PATH"

