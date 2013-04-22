#!/usr/bin/env nextflow

channel = new Channel(1,2,3)

task {
   input channel
   threads 4
   errorStrategy 'ignore'

   "echo $channel; exit 1"

 }

 sleep 500

 channel2 = new Channel(4,5,6)

 task {
   input channel2
   threads 4

   "echo $channel; exit 1"

 }