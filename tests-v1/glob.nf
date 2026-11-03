#!/usr/bin/env nextflow
/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

process recurseDir {

    output:
    file 'folder/**.fa'
    file 'folder/**/*.txt' 

    """
    mkdir -p folder/x
    mkdir -p folder/y
    touch folder/file1.txt
    touch folder/file2.fa
    touch folder/x/file3.fa
    touch folder/y/file4.fa
    touch folder/y/file5.txt
    """
}

workflow {
  recurseDir()
  recurseDir.out[0] | flatten | view { "result1: " + it.name }
  recurseDir.out[1] | flatten | view { "result2: " + it.name }
}
