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

package misc
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Session
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


def session = new Session()

def values = new DataflowQueue()
def strings = new DataflowVariable()

values << 1 << 2 << PoisonPill.instance
strings << 'a'


def (stdout, aln) = session.createProcessor()
        .name('step1')
        .input( x: values, y: strings )
        .output('-')
        .output('file')
        .script {
            """\
            echo 'Hello world! - ${x}'
            echo "some content $y" > file
            """
        }
        .run()



session.createProcessor()
        .name('step2')
        .input(str: aln)
        .script {
                """
                echo $str
                """
            }
        .run()


session.createProcessor()
        .name('step3')
        .input(file: aln)
        .script{
                """
                cat $fInputFile
                """
            }
        .run()



session.terminate()
