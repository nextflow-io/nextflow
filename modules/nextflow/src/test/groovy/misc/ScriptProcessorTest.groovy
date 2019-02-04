/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
