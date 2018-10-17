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
 
seqs = Channel.fromPath("$baseDir/data/p{1,2,3}.fa")

process foo {
    echo true

    input:
    file 'dir1/link_*.fasta' from seqs.toList()

    output:
    file 'dir2/*' into result mode flatten

    '''
    ls dir1 | sort
    mkdir dir2
    echo hello > dir2/alpha.txt
    echo world > dir2/beta.txt
    '''
}

result.println { it.name }
