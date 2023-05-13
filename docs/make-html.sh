#!/bin/bash

docker run -v $(pwd):/tmp nextflow/sphinx:5.3.0 -- make html
echo "Done. See _build/html/index.html"
