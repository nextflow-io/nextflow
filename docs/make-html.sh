#!/bin/bash

docker run -v $(pwd):/tmp nextflow/sphinx:5.3.0 -- make html
