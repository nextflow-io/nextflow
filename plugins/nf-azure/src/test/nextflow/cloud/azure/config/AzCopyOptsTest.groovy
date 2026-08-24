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

package nextflow.cloud.azure.config

import spock.lang.Specification

/**
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
class AzCopyOptsTest extends Specification {

    def 'should get block size and blob tier'() {
        when:
        def opts1 = new AzCopyOpts([:])
        then:
        opts1.blobTier == 'None'
        opts1.blockSize == '4'

        when:
        def opts2 = new AzCopyOpts(
                [blobTier: 'Hot', blockSize: '100'])
        then:
        opts2.blobTier == 'Hot'
        opts2.blockSize == '100'

    }

}
