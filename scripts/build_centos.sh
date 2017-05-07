#!/bin/bash
set -e

#based on
#https://sourceforge.net/p/qwwad/wiki/Installing%20from%20a%20package/

# This script builds qwwad on the centos 6 platform
# This successfully builds in a docker image with all of the
# following dependencies installed
# To be run from the project root

yum install -y epel-release
yum update -y
yum install -y gcc gcc-c++ gcc-gfortran help2man armadillo-devel boost148-program-options boost148-devel gsl-devel lapack-devel libxml++-devel cmake3 git wget

if [ ! -f /etc/yum.repos.d/devtools-2.repo ]; then \
  wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
fi
yum install -y devtoolset-2-gcc devtoolset-2-binutils devtoolset-2-gcc-c++ devtoolset-2-gcc-gfortran

#scl enable devtoolset-2 bash
export CMAKE_CXX_COMPILER=/opt/rh/devtoolset-2/root/usr/bin/g++
export CMAKE_FORTRAN_COMPILER=/opt/rh/devtoolset-2/root/usr/bin/gfortran
export CMAKE_C_COMPILER=/opt/rh/devtoolset-2/root/usr/bin/gcc
export CMAKE_INSTALL_PREFIX=${HOME}/qwwad



#cd ${HOME}
#git clone https://github.com/QWWAD/qwwad.git
#cd qwwad
ln -sf /usr/include/boost148/boost /usr/include/boost
ln -sf /usr/lib64/libboost_program_options-mt.so.1.48.0 /usr/lib64/libboost_program_options-mt.so
ln -sf /usr/lib64/libboost_program_options.so.1.48.0 /usr/lib64/libboost_program_options.so

(rm -rf build && mkdir build && cd build && cmake3 -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ..)
(cd build && make -j2)
(cd build && make install)

