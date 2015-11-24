#!/bin/bash

export NXF_CMD=$PWD/nextflow; 
(
 $NXF_CMD info
 cd validation
 bash -x test.sh
)
