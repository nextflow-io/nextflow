#!/bin/bash

docker run -v $(pwd):/tmp $(wave -f Dockerfile --context .) -- make html
echo "Done. See _build/html/index.html"
