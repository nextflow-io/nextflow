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

package nextflow.cloud.aws.batch

import java.nio.file.Paths

import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsOptionsTest extends Specification {

    def 'should return aws cli' () {

        given:
        AwsOptions opts

        when:
        opts = new AwsOptions()
        then:
        opts.awsCli == 'aws'

        when:
        opts = new AwsOptions(cliPath: '/foo/bin/aws')
        then:
        opts.awsCli == '/foo/bin/aws'

        when:
        opts = new AwsOptions(cliPath: '/foo/bin/aws', region: 'eu-west-1')
        then:
        opts.awsCli == '/foo/bin/aws --region eu-west-1'
    }

    def 'should get max connection'  () {
        given:
        def sess = Mock(Session)  {
            getConfig() >> [aws:[batch:[maxParallelTransfers: 5]]]
        }
        AwsOptions opts

        when:
        opts = new AwsOptions()
        then:
        opts.maxParallelTransfers == AwsOptions.MAX_TRANSFER

        when:
        opts = new AwsOptions(sess)
        then:
        opts.maxParallelTransfers == 5

    }

    def 'should get aws options' () {
        given:
        def sess = Mock(Session)  {
            getConfig() >> [aws:
                                [
                                    batch:[
                                        cliPath: '/foo/bar/aws',
                                        maxParallelTransfers: 5,
                                        maxTransferAttempts: 3,
                                        delayBetweenAttempts: '9 sec',
                                        jobRole: 'aws::foo::bar',
                                        volumes: '/foo,/this:/that'],
                                    client: [
                                        uploadStorageClass: 'my-store-class',
                                        storageEncryption: 'my-ecrypt-class'],
                                    region: 'aws-west-2'
                                ]
            ]
        }

        def exec = Mock(AwsBatchExecutor)
        exec.getSession() >> sess
        exec.getRemoteBinDir() >> Paths.get('/remote/bin/path')

        when:
        def opts = new AwsOptions(sess)
        then:
        opts.maxParallelTransfers == 5
        opts.maxTransferAttempts == 3 
        opts.delayBetweenAttempts.seconds == 9
        opts.storageClass == 'my-store-class'
        opts.storageEncryption == 'my-ecrypt-class'
        opts.region == 'aws-west-2'
        opts.jobRole == 'aws::foo::bar'
        opts.volumes == ['/foo','/this:/that']

        when:
        opts = new AwsOptions(exec)
        then:
        opts.remoteBinDir == '/remote/bin/path'

    }

    def 'should parse volumes list' () {

        given:
        def executor = Spy(AwsOptions)

        expect:
        executor.makeVols(OBJ) == EXPECTED

        where:
        OBJ             | EXPECTED
        null            | []
        'foo'           | ['foo']
        'foo, bar'      | ['foo','bar']
        '/foo/,/bar///' | ['/foo','/bar']
        ['/this','/that'] | ['/this','/that']
        ['/foo/bar/']   | ['/foo/bar']

    }


    @Unroll
    def 'should return aws options'() {
        given:
        def cfg = [
                aws: [client: [
                        uploadStorageClass: awsStorClass,
                        storageEncryption  : awsStorEncrypt]],
                executor: [
                        awscli: awscliPath
                ]
        ]
        def session = new Session(cfg)

        when:
        def opts = new AwsOptions(session)
        then:
        opts.cliPath == awscliPath
        opts.storageClass == awsStorClass
        opts.storageEncryption == awsStorEncrypt

        where:
        awscliPath      | awsStorClass | awsStorEncrypt
        null            | null         | null
        '/foo/bin/aws'  | 'STANDARD'   | 'AES256'

    }

    def 'should validate aws options' () {

        when:
        def opts = new AwsOptions()
        then:
        opts.getCliPath() == null
        opts.getStorageClass() == null
        opts.getStorageEncryption() == null

        when:
        opts = new AwsOptions(cliPath: '/foo/bin/aws', storageClass: 'STANDARD', storageEncryption: 'AES256')
        then:
        opts.getCliPath() == '/foo/bin/aws'
        opts.getStorageClass() == 'STANDARD'
        opts.getStorageEncryption() == 'AES256'

        when:
        opts = new AwsOptions(storageClass: 'foo')
        then:
        opts.getStorageClass() == null

        when:
        opts = new AwsOptions(storageEncryption: 'abr')
        then:
        opts.getStorageEncryption() == null

        when:
        new AwsOptions(cliPath: 'bin/aws')
        then:
        thrown(ProcessUnrecoverableException)

        when:
        new AwsOptions(cliPath: '/foo/aws')
        then:
        thrown(ProcessUnrecoverableException)
    }

    def 'should add a volume' () {
        given:
        def opts = new AwsOptions()

        when:
        opts.addVolume(Paths.get('/some/dir'))
        then:
        opts.volumes == ['/some/dir']

        when:
        opts.addVolume(Paths.get('/other/dir'))
        opts.addVolume(Paths.get('/other/dir'))
        then:
        opts.volumes == ['/some/dir', '/other/dir']
    }

}
