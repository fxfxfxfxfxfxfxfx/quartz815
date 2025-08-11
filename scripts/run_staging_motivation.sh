#!/bin/bash

QUARTZ_HOME=$(pwd)
export PATH=$QUARTZ_HOME/external/HiGHS/build/bin/:$PATH

./build/benchmark_snuqs_num_stages
