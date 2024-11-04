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

package nextflow.executor


import com.google.api.client.http.HttpResponseException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderWithGoogleTest extends Specification {

    def 'should check retryable errors' () {
        expect:
        BashWrapperBuilder.isRetryable0(ERROR) == EXPECTED
        where:
        ERROR                                                                   | EXPECTED
        new HttpResponseException(GroovyMock(HttpResponseException.Builder))    | true
        new Exception()                                                         | false
    }

}
