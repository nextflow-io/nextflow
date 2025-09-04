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

// Test multi-process entry selection with parameter mapping
process preprocessData {
    debug true
    
    input:
    path inputFile
    val quality
    
    script:
    """
    echo "Multi-process test: preprocessing"
    echo "Input file: ${inputFile}"
    echo "Quality threshold: ${quality}"
    
    if [ -f "${inputFile}" ]; then
        echo "Preprocessing ${inputFile} with quality ${quality}"
        head -n 5 "${inputFile}"
    else
        echo "Input file ${inputFile} not found"
    fi
    
    echo "Preprocessing completed"
    """
}

process analyzeResults {
    debug true
    
    input:
    val experimentId
    path resultsFile
    val mode
    
    script:
    """
    echo "Multi-process test: analysis"
    echo "Experiment ID: ${experimentId}"
    echo "Results file: ${resultsFile}"
    echo "Analysis mode: ${mode}"
    
    if [ -f "${resultsFile}" ]; then
        echo "Analyzing results for experiment ${experimentId} in ${mode} mode"
        tail -n 3 "${resultsFile}"
    else
        echo "Results file ${resultsFile} not found"
    fi
    
    echo "Analysis completed for experiment ${experimentId}"
    """
}

process generateReport {
    debug true
    
    input:
    val reportTitle
    path dataPath
    
    script:
    """
    echo "Multi-process test: reporting"
    echo "Report title: ${reportTitle}"
    echo "Data path: ${dataPath}"
    
    if [ -f "${dataPath}" ]; then
        echo "Generating report '${reportTitle}' from ${dataPath}"
        echo "Data file size:"
        wc -c "${dataPath}"
    else
        echo "Data path ${dataPath} not found"
    fi
    
    echo "Report generation completed"
    """
}