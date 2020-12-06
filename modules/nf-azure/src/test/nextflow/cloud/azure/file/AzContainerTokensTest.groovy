/*
 * Copyright 2020, Microsoft Corp
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

package nextflow.cloud.azure.file

import nextflow.cloud.azure.file.AzContainerTokens
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzContainerTokensTest extends Specification {

    @Unroll
    def 'should parse az uri' () {
        expect:
        AzContainerTokens.parse(AZ_URI) == TOKENS

        where:
        AZ_URI                                              | TOKENS
        'azb://foo:/'                                       | new AzContainerTokens('foo:', '/', null)
        'azb://foo:/this/that?account=nfstore'              | new AzContainerTokens('foo:', '/this/that', 'nfstore')
        'azb://my-data:/nf-DE8YJPrG.txt?account=nfstore'    | new AzContainerTokens('my-data:', '/nf-DE8YJPrG.txt', 'nfstore')
    }
}
