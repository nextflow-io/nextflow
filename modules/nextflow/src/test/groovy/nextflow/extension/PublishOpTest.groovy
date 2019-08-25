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

package nextflow.extension

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.TimeoutException

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import test.BaseSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PublishOpTest extends BaseSpec {

    def 'should infer task dir' () {
        given:
        def BASE = '/some/work/dir' as Path
        def BUCKET = 's3:/other/bucket/dir' as Path
        def sess = Mock(Session) {
            getWorkDir() >> BASE
            getBucketDir() >> BUCKET
        }
        Global.session = sess

        def op = new PublishOp(Mock(DataflowReadChannel), [to:'/target'])

        when:
        def result = op.getTaskDir( BASE.resolve('xx/yyyy/this/and/that.txt') )
        then:
        result == Paths.get('/some/work/dir/xx/yyyy')

        when:
        result = op.getTaskDir( BUCKET.resolve('pp/qqqq/other/file.fasta') )
        then:
        result == 's3:/other/bucket/dir/pp/qqqq' as Path


        when:
        result = op.getTaskDir( BASE.resolve('xx/foo.txt') )
        then:
        result == null
    }

    def 'should publish files' () {
        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('file1.txt'); file1.text = 'Hello'
        def file2 = folder.resolve('file2.txt'); file2.text = 'world'
        def target = folder.resolve('target/dir')


        def BASE = folder
        def sess = Mock(Session) {
            getWorkDir() >> BASE
        }
        Global.session = sess

        and:
        def ch = new DataflowQueue()
        ch.bind(file1)
        ch.bind(file2)
        ch.bind(Channel.STOP)

        when:
        def now = System.currentTimeMillis()
        def op = new PublishOp(ch, [to:target, mode:'symlink']) .apply()
        while( !op.complete ) { sleep 100; if( System.currentTimeMillis()-now>5_000) throw new TimeoutException() }
        then:
        target.resolve('file1.txt').text == 'Hello'
        target.resolve('file2.txt').text == 'world'

        cleanup:
        folder?.deleteDir()
    }

}
