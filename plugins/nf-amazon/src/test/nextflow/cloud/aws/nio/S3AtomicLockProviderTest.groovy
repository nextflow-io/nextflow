/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.aws.nio

import java.nio.file.Path
import java.nio.file.ProviderMismatchException

import spock.lang.Specification

class S3AtomicLockProviderTest extends Specification {

    def 'handles only the s3 scheme'() {
        given:
        def provider = new S3AtomicLockProvider()
        expect:
        provider.canHandle('s3')
        !provider.canHandle('gs')
        !provider.canHandle('az')
    }

    def 'tryCreate delegates to the conditional PUT'() {
        given:
        def client = Mock(S3Client)
        def fs = Mock(S3FileSystem) { getClient() >> client }
        def path = Mock(S3Path) { getFileSystem() >> fs; getBucket() >> 'b'; getKey() >> 'k' }
        def provider = new S3AtomicLockProvider()

        when:
        def first = provider.tryCreate(path)
        then:
        1 * client.putObjectIfAbsent('b','k') >> true
        first

        when:
        def second = provider.tryCreate(path)
        then:
        1 * client.putObjectIfAbsent('b','k') >> false
        !second
    }

    def 'tryCreate throws for a path of another provider, instead of reporting a lost race'() {
        when:
        new S3AtomicLockProvider().tryCreate(Mock(Path))
        then: 'a false return would be read as "lost the race" and retried forever'
        thrown(ProviderMismatchException)
    }

}
