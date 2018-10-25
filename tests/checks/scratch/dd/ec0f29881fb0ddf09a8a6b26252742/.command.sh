#!/bin/bash -ue
awk '/^>/{f="seq_"++d} {print > f}' < input.fa
