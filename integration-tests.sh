#!/bin/bash

export NXF_CMD=$PWD/nextflow; 
(
 cd validation
 bash -x test.sh
)
