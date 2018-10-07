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

/*
 * Show how get data from the 'stdin' virtual file
 *
 * The data is supposed to be a FASTA file which is splitted
 * is chunks containing a single sequence
 *
 * The number of sequences in each each is controlled by the command
 * line parameter '--chunkSize' (--chunk-size is a synonym for the same)
 */


params.chunkSize = 1
sequences = Channel.create()

stdin.splitFasta( by: params.chunkSize, into: sequences )

process foo {
    echo true

    input:
    stdin sequences

    "cat -"
}
