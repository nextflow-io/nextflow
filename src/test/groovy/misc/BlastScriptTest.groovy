package misc
/*
 * Copyright (c) 2012, the authors.
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

import nextflow.Nextflow
import nextflow.Session
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


def session = new Session()

def fasta='/Users/ptommaso/sample.fa'
def DB='/Users/ptommaso/tools/blast-db/pdb/pdb'

def sequences = session.createProcessor()
        .input(fasta:[fasta])
        .output('seq_*')
        .script {
            """
            gcsplit $fasta '%^>%' '/^>/' '{*}' -f seq_
            """
        }
        .run()


def singleBlast = session.createProcessor()
    .input(seq: sequences)
    .output('blast_result')
    .script {
        """
        blastp -db $DB -query $seq -outfmt 6 > blast_result
        """
    }
    .run()


def allBlast = session.createProcessor()
    .input(blast:singleBlast)
    .output('allBlast')
    .setSharedWorkDirectory(true)
    .setBindOnTermination(true)
    .script {
        """
        cat $blast >> allBlast
        """
     }
    .run()

def result = session.createProcessor()
    .input(result:allBlast)
    .script {
        """
        sort $result
        """
        }
    .run()


println Nextflow.read(result)

session.terminate()