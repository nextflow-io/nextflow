#!/bin/bash -ue
gcsplit query.fa '%^>%' '/^>/' '{*}' -f seq_
