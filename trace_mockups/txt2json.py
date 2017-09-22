#!/usr/bin/env python
# Convert Nextflow text file to JSON format

import json

headers = []
data = {'trace': []}
with open ('NGI-RNAseq_trace.txt') as fh:
    for l in fh:
        s = l.strip().split("\t")
        if len(headers) == 0:
            headers = s
        else:
            d = {}
            for i, v in enumerate(s):
                try:
                    d[headers[i]] = float(v)
                except ValueError:
                    d[headers[i]] = v
            data['trace'].append(d)

print(json.dumps(data, indent=4))
