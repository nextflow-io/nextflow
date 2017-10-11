#!/usr/bin/env bash

VER='16.0.1'
FILE="guava-$VER.jar"
TARGET="guava-cache-$VER.jar"

if [ ! -f $FILE ]; then 
wget -O $FILE http://search.maven.org/remotecontent?filepath=com/google/guava/guava/16.0.1/$FILE
fi

java -jar jarjar-1.4.jar process guava.rules $FILE ../lib/$TARGET
