#!/bin/bash

export NXF_CMD=$PWD/launch.sh; 
(
 cd validation
 bash -x test.sh
)
