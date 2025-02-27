/*
 * Copyright 2013-2024, Seqera Labs
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

import com.amazonaws.services.s3.model.CannedAccessControlList
import nextflow.Session
import nextflow.cloud.aws.config.AwsConfig
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
        opts = new AwsOptions(awsConfig: new AwsConfig([:]))
        then:
        opts.awsCli == 'aws'

        when:
        opts = new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/foo/bin/aws']))
        then:
        opts.awsCli == '/foo/bin/aws'

        when:
        opts = new AwsOptions(awsConfig: new AwsConfig(region: 'eu-west-1', batch: [cliPath: '/foo/bin/aws']))
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
        opts = new AwsOptions(awsConfig: new AwsConfig([:]))
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
                                        cliPath: '/foo/bin/aws',
                                        maxParallelTransfers: 5,
                                        maxTransferAttempts: 3,
                                        delayBetweenAttempts: '9 sec',
                                        jobRole: 'aws::foo::bar',
                                        volumes: '/foo,/this:/that'],
                                    client: [
                                        uploadStorageClass: 'STANDARD',
                                        storageEncryption: 'AES256'],
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
        opts.storageClass == 'STANDARD'
        opts.storageEncryption == 'AES256'
        opts.region == 'aws-west-2'
        opts.jobRole == 'aws::foo::bar'
        opts.volumes == ['/foo','/this:/that']

        when:
        opts = new AwsOptions(exec)
        then:
        opts.remoteBinDir == '/remote/bin/path'

    }

    def 'should set aws kms key' () {
        when:
        def sess1 = Mock(Session)  {
            getConfig() >> [aws: [ client: [ storageKmsKeyId: 'my-kms-key']]]
        }
        and:
        def opts = new AwsOptions(sess1)
        then:
        opts.storageKmsKeyId == 'my-kms-key'
        opts.storageEncryption == null

        when:
        def sess2 = Mock(Session)  {
            getConfig() >> [aws: [ client: [ storageKmsKeyId: 'my-kms-key', storageEncryption: 'aws:kms']]]
        }
        and:
        def opts2 = new AwsOptions(sess2)
        then:
        opts2.storageKmsKeyId == 'my-kms-key'
        opts2.storageEncryption == 'aws:kms'    // <-- allow explicit `storageEncryption`

    }



    @Unroll
    def 'should return aws options'() {
        given:
        def cfg = [
                aws: [client: [
                        uploadStorageClass: awsStorClass,
                        storageEncryption  : awsStorEncrypt],
                      batch: [ cliPath: awscliPath ]]
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
        def opts = new AwsOptions(awsConfig: new AwsConfig([:]))
        then:
        opts.getCliPath() == null
        opts.getStorageClass() == null
        opts.getStorageEncryption() == null

        when:
        opts = new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/foo/bin/aws'], client: [storageClass: 'STANDARD', storageEncryption: 'AES256']))
        then:
        opts.getCliPath() == '/foo/bin/aws'
        opts.getStorageClass() == 'STANDARD'
        opts.getStorageEncryption() == 'AES256'

        when:
        opts = new AwsOptions(awsConfig: new AwsConfig(client:[storageClass: 'foo']))
        then:
        opts.getStorageClass() == null

        when:
        opts = new AwsOptions(awsConfig: new AwsConfig(client:[storageEncryption: 'abr']))
        then:
        opts.getStorageEncryption() == null

        when:
        opts = new AwsOptions(awsConfig: new AwsConfig(client:[storageKmsKeyId: 'arn:aws:kms:eu-west-1:1234567890:key/e97ecf28-951e-4700-bf22-1bd416ec519f']))
        then:
        opts.getStorageKmsKeyId() == 'arn:aws:kms:eu-west-1:1234567890:key/e97ecf28-951e-4700-bf22-1bd416ec519f'

        when:
        new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: 'bin/aws']))
        then:
        thrown(ProcessUnrecoverableException)

        when:
        new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/foo/aws']))
        then:
        thrown(ProcessUnrecoverableException)
    }

    def 'should add a volume' () {
        given:
        def opts = new AwsOptions(awsConfig: new AwsConfig([:]))

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

    @Unroll
    def 'should get aws cli path' () {
        def session = new Session(CONFIG)

        when:
        def opts = new AwsOptions(session)
        then:
        opts.cliPath == S3CLI_PATH
        opts.s5cmdPath == S5CMD_PATH

        where:
        CONFIG                                                                      | S3CLI_PATH        | S5CMD_PATH
        [aws:[batch:[:]]]                                                           | null              | null
        [aws:[batch:[cliPath: '/usr/bin/aws']]]                                     | '/usr/bin/aws'    | null
        [aws:[batch:[cliPath: 's5cmd']]]                                            | null              | null
        [aws:[batch:[platformType: 'fargate', cliPath: 's5cmd']]]                   | null              | 's5cmd'
        [aws:[batch:[platformType: 'fargate', cliPath: '/some/path/s5cmd']]]        | null              | '/some/path/s5cmd'
        [aws:[batch:[platformType: 'fargate', cliPath: 's5cmd --foo']]]             | null              | 's5cmd --foo'
        [aws:[batch:[platformType: 'fargate', cliPath: '/some/path/s5cmd --foo']]]  | null              | '/some/path/s5cmd --foo'
    }
    
    def 'should parse s3 acl' ( ) {
        when:
        def opts = new AwsOptions(new Session(aws:[client:[s3Acl: 'PublicRead']]))
        then:
        opts.getS3Acl() == CannedAccessControlList.PublicRead


        when:
        opts = new AwsOptions(new Session(aws:[client:[s3Acl: 'public-read']]))
        then:
        opts.getS3Acl() == CannedAccessControlList.PublicRead


        when:
        opts = new AwsOptions(new Session(aws:[client:[s3Acl: 'unknown']]))
        then:
        thrown(IllegalArgumentException)
    }

}
