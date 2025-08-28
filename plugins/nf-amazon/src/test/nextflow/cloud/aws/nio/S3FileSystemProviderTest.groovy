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
 *
 */

package nextflow.cloud.aws.nio

import software.amazon.awssdk.services.s3.model.ObjectCannedACL
import software.amazon.awssdk.services.s3.model.ServerSideEncryption
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3FileSystemProviderTest extends Specification {

    def 'should create filesystem from config'(){
        given:
        def config = [
            client: [
                anonymous: true,
                s3Acl: 'Private',
                connectionTimeout: 20000,
                endpoint: 'https://s3.eu-west-1.amazonaws.com',
                minimumPartSize: '7MB',
                multipartThreshold: '32MB',
                maxConnections: 100,
                maxErrorRetry: 3,
                socketTimeout: 20000,
                requesterPays: true,
                s3PathStyleAccess: true,
                proxyHost: 'host.com',
                proxyPort: 80,
                proxyScheme: 'https',
                proxyUsername: 'user',
                proxyPassword: 'pass',
                storageEncryption: 'AES256',
                storageKmsKeyId: 'arn:key:id',
                uploadMaxThreads: 15,
                uploadChunkSize: '7MB',
                uploadMaxAttempts: 4,
                uploadRetrySleep: '200ms'
            ],
            accessKey: '123456abc',
            secretKey: '78910def',
            profile: 'test'
        ]
        def provider = new S3FileSystemProvider();
        when:
        def fs = provider.newFileSystem(new URI("s3:///bucket/key"), config) as S3FileSystem
        then:
        fs.getBucketName() == 'bucket'
        def client = fs.getClient()
        client.client != null
        client.cannedAcl == ObjectCannedACL.PRIVATE
        client.storageEncryption == ServerSideEncryption.AES256
        client.isRequesterPaysEnabled == true
        client.kmsKeyId == 'arn:key:id'
        client.factory.accessKey() == '123456abc'
        client.factory.secretKey() == '78910def'
        client.factory.profile() == 'test'
        client.factory.config.s3Config.anonymous == true
        client.factory.config.s3Config.endpoint == 'https://s3.eu-west-1.amazonaws.com'
        client.factory.config.s3Config.pathStyleAccess == true
        fs.properties().getProperty('proxy_host') == 'host.com'
        fs.properties().getProperty('proxy_port') == '80'
        fs.properties().getProperty('proxy_scheme') == 'https'
        fs.properties().getProperty('proxy_username') == 'user'
        fs.properties().getProperty('proxy_password') == 'pass'
        fs.properties().getProperty('socket_timeout') == '20000'
        fs.properties().getProperty('connection_timeout') == '20000'
        fs.properties().getProperty('max_connections') == '100'
        fs.properties().getProperty('max_error_retry') == '3'
        fs.properties().getProperty('upload_max_attempts') == '4'
        fs.properties().getProperty('upload_retry_sleep') == '200'
        fs.properties().getProperty('upload_chunk_size') == '7340032' //7MB
        fs.properties().getProperty('upload_max_threads') == '15'
        fs.properties().getProperty('minimum_part_size') == '7340032' //7MB
        fs.properties().getProperty('multipart_threshold') == '33554432' //32MB
    }


}
