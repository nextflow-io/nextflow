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

process recurseDir {

    output:
    file 'folder/**.fa' into result1
    file 'folder/**/*.txt' into result2

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


result1.flatten().subscribe { println "result1: " + it.name }
result2.flatten().subscribe { println "result2: " + it.name }
