/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.script

import spock.lang.Specification

/**
 * Tests for {@link PlatformMetadata}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PlatformMetadataTest extends Specification {

    def 'should create with default constructor'() {
        when:
        def meta = new PlatformMetadata()

        then:
        meta.workflowId == null
    }

    def 'should create with workflowId'() {
        when:
        def meta = new PlatformMetadata('abc123')

        then:
        meta.workflowId == 'abc123'
    }

    def 'should allow setting workflowId after construction'() {
        given:
        def meta = new PlatformMetadata()

        when:
        meta.workflowId = 'xyz789'

        then:
        meta.workflowId == 'xyz789'
    }

    def 'should create with null workflowUrl by default'() {
        when:
        def meta = new PlatformMetadata()

        then:
        meta.workflowUrl == null
    }

    def 'should allow setting workflowUrl after construction'() {
        given:
        def meta = new PlatformMetadata()

        when:
        meta.workflowUrl = 'https://cloud.seqera.io/watch/abc123'

        then:
        meta.workflowUrl == 'https://cloud.seqera.io/watch/abc123'
    }
}
