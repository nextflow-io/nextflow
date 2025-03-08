#!/bin/bash
launch_cmd=../../launch.sh
export NXF_ANSI_LOG=false
export WITH_DOCKER=${WITH_DOCKER:=''}
export NXF_CMD=$(realpath $launch_cmd)
export TEST_JDK=17
bash run.sh "$1"
