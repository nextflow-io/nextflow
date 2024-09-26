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
package nextflow.util

import spock.lang.Specification
/**
 * Common task utilities.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class TaskUtilsTest extends Specification{

    def 'should include the tower prefix'() {
        given:
        def name = 'job_1'

        expect:
        TaskUtils.includeTowerPrefix(name, ENV) == EXPECTED

        where:
        ENV                         | EXPECTED
        [TOWER_WORKFLOW_ID: '1234'] | "tw-1234-job_1"
        [:]                         | "job_1"
    }
}
