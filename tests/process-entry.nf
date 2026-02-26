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

// Test single process auto-execution with parameter mapping
process analyzeData {
    debug true
    
    input:
    val sampleId
    path dataFile
    val threads
    
    script:
    """
    echo "Single process execution test"
    echo "Sample ID: ${sampleId}"
    echo "Data file: ${dataFile}"
    echo "Threads: ${threads}"
    
    if [ -f "${dataFile}" ]; then
        echo "Processing file with ${threads} threads..."
        wc -l "${dataFile}"
    else
        echo "Warning: File ${dataFile} not found"
    fi
    
    echo "Analysis completed for sample ${sampleId}"
    """
}