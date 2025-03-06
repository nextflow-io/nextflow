launch_cmd=../../launch.sh
export NXF_ANSI_LOG=false
export WITH_DOCKER=-with-docker
export NXF_CMD=$(realpath $launch_cmd)
export NXF_IGNORE_WARN_DSL2=true
export TEST_JDK=11
bash run.sh "$1"
