/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow

import java.nio.file.Files

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GlobalTest extends Specification {


    def testAwsCredentials() {

        expect:
        Global.getAwsCredentials0(null, null) == null
        Global.getAwsCredentials0([AWS_ACCESS_KEY: 'x', AWS_SECRET_KEY: '222'], null) == ['x','222']
        Global.getAwsCredentials0([AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999'], null) == ['q','999']
        Global.getAwsCredentials0([AWS_ACCESS_KEY: 'x', AWS_SECRET_KEY: '222',  AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999'], null) == ['q','999']

        Global.getAwsCredentials0([AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999'], [aws:[accessKey: 'b', secretKey: '333']]) == ['b','333']
        Global.getAwsCredentials0(null, [aws:[accessKey: 'b', secretKey: '333']]) == ['b','333']
        Global.getAwsCredentials0(null, [aws:[accessKey: 'b']]) == null

    }

    def testAwsCredentialsWithFile() {

        given:
        def file = Files.createTempFile('test','test')
        file.text = '''
            [default]
            aws_access_key_id = aaa
            aws_secret_access_key = bbbb
            '''

        Global.getAwsCredentials0(null, null, [file]) == ['aaa','bbbb']
        Global.getAwsCredentials0([AWS_ACCESS_KEY: 'x', AWS_SECRET_KEY: '222'], null, [file]) == ['x','222']

        cleanup:
        file?.delete()

    }

    def 'should normalize aws config' () {

        given:
        def config = [uploadMaxThreads: 5, uploadChunkSize: 1000, uploadStorageClass: 'STANDARD' ]
        when:
        def norm = Global.normalizeAwsClientConfig(config)
        then:
        norm.upload_storage_class == 'STANDARD'
        norm.upload_chunk_size == '1000'
        norm.upload_max_threads == '5'

        when:
        config.uploadChunkSize = '10MB'
        then:
        Global.normalizeAwsClientConfig(config).upload_chunk_size == '10485760'

        when:
        config.uploadChunkSize = '1024'
        then:
        Global.normalizeAwsClientConfig(config).upload_chunk_size == '1024'

        when:
        config.uploadChunkSize = new MemoryUnit('2 MB')
        then:
        Global.normalizeAwsClientConfig(config).upload_chunk_size == '2097152'

        when:
        config.uploadRetrySleep = '10 sec'
        then:
        Global.normalizeAwsClientConfig(config).upload_retry_sleep == '10000'

        when:
        config.uploadRetrySleep = Duration.of('5 sec')
        then:
        Global.normalizeAwsClientConfig(config).upload_retry_sleep == '5000'
    }


}
