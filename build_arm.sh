#!/bin/bash

curr_path=$(pwd)
echo "curr_path=$curr_path"

# COMPILER
export PATH=${PATH}:/opt/toolchain/7.5.0/gcc-linaro-7.5.0-2019.12-x86_64_aarch64-linux-gnu/bin
GCC_COMPILER=aarch64-linux-gnu

# build
BUILD_DIR=./build
if [[ ! -d "${BUILD_DIR}" ]]; then
  mkdir -p ${BUILD_DIR}
fi
cd ${BUILD_DIR}


cmake .. \
    -DCMAKE_C_COMPILER=${GCC_COMPILER}-gcc \
    -DCMAKE_CXX_COMPILER=${GCC_COMPILER}-g++ \
    -DCMAKE_INSTALL_PREFIX=../install \
    -DCMAKE_BUILD_TYPE=Release \
    
make -j8
make install
