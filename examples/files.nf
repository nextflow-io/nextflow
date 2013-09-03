/*
 * Copyright (c) 2013, the authors.
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



params.in = '~/sample.fa'
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

task {
    input file:'query.fa', from: file(params.in)
    output file: 'seq_*', into: splits, joint: true

    """
    $SPLIT query.fa '%^>%' '/^>/' '{*}' -f seq_
    """
}


task{
    input file: 'chunk', from: splits
    echo true

    """
    cat chunk* | rev
    """
}
