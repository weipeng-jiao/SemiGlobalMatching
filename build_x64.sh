#!/bin/bash

cur_dir=$(cd `dirname $0`; pwd)
echo cur_dir=$cur_dir
BUILD_DIR=${cur_dir}/build
if [[ ! -d "${BUILD_DIR}" ]]; then
  mkdir -p ${BUILD_DIR}
fi
cd ${BUILD_DIR}
cmake .. \
    -DCMAKE_INSTALL_PREFIX=../install   
make -j8
make install
