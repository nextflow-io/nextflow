#!/usr/bin/env nextflow 
/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

params.in = "$baseDir/data/sample.fa"

/*
 * Splits the input file in chunks containing a single sequences,
 * and send each of it over the 'seq' channel
 */
seq = Channel.fromPath(params.in).splitFasta()

/*
 * For each sequence that is sent over the 'seq' channel
 * the below task is executed
 */
process ampaTask {

    input:
    file seq

    output:
    file result

    // The BASH script to be executed - for each - sequence
    """
    AMPA.pl -in=${seq} -noplot -rf=result -df=data
    """

}

/*
 * print out each 'result' produced by the above step
 */
result.subscribe { println it.text }
