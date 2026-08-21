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

package nextflow.processor

import java.nio.file.NoSuchFileException

import nextflow.util.RetryConfig
import spock.lang.Specification
class TaskFileRetryTest extends Specification {

    private static final RetryConfig FAST_RETRY = new RetryConfig([delay: '20ms', maxDelay: '50ms', maxAttempts: 5])

    def 'should return the result on the first successful attempt' () {
        when:
        def result = TaskFileRetry.withRetry(FAST_RETRY, [NoSuchFileException], 'test file') {
            'hello'
        }
        then:
        result == 'hello'
    }

    def 'should retry until the body succeeds' () {
        given:
        int attempts = 0

        when:
        def result = TaskFileRetry.withRetry(FAST_RETRY, [NoSuchFileException], 'test file') {
            if( ++attempts < 3 )
                throw new NoSuchFileException('missing')
            return 'hello'
        }
        then:
        result == 'hello'
        attempts == 3
    }

    def 'should give up after maxAttempts and rethrow the original exception' () {
        given:
        int attempts = 0

        when:
        TaskFileRetry.withRetry(FAST_RETRY, [NoSuchFileException], 'test file') {
            ++attempts
            throw new NoSuchFileException('missing')
        }
        then:
        def e = thrown(NoSuchFileException)
        e.message == 'missing'
        attempts == FAST_RETRY.maxAttempts
    }

    def 'should not retry an exception outside the retryOn list' () {
        given:
        int attempts = 0

        when:
        TaskFileRetry.withRetry(FAST_RETRY, [NoSuchFileException], 'test file') {
            ++attempts
            throw new IllegalStateException('boom')
        }
        then:
        def e = thrown(IllegalStateException)
        e.message == 'boom'
        attempts == 1
    }

}
