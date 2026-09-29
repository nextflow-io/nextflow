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

package io.seqera.config

import spock.lang.Specification

/**
 * Unit tests for SchedulingRequirementOpts
 *
 * @author Jon Martí <jonathan.marti@seqera.io>
 */
class SchedulingRequirementOptsTest extends Specification {

    def 'should create with empty config' () {
        when:
        def opts = new SchedulingRequirementOpts([:])

        then:
        opts.maxCpusPerUser == null
    }

    def 'should create with maxCpusPerUser' () {
        when:
        def opts = new SchedulingRequirementOpts([maxCpusPerUser: 16])

        then:
        opts.maxCpusPerUser == 16
    }

    def 'should coerce a string maxCpusPerUser to integer' () {
        when:
        def opts = new SchedulingRequirementOpts([maxCpusPerUser: '16'])

        then:
        opts.maxCpusPerUser == 16
    }

    def 'should reject a non-positive maxCpusPerUser' () {
        when:
        new SchedulingRequirementOpts([maxCpusPerUser: value])

        then:
        def err = thrown(IllegalArgumentException)
        err.message.contains("'seqera.executor.schedulingRequirement.maxCpusPerUser'")
        err.message.contains('positive integer')

        where:
        value << [0, -1, -16]
    }
}
