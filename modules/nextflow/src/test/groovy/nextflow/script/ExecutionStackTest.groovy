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

package nextflow.script

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutionStackTest extends Specification {

    def setupSpec() {
        ExecutionStack.reset()
    }

    def 'should verify push and pop semantics' () {

        given:
        def s1 = Mock(BaseScript)
        def s2 = Mock(WorkflowDef)
        def s3 = Mock(WorkflowDef)

        expect:
        ExecutionStack.size()==0

        when:
        ExecutionStack.push(s1)
        then:
        ExecutionStack.size()==1
        ExecutionStack.current() == s1
        !ExecutionStack.withinWorkflow()
        when:
        ExecutionStack.push(s2)
        then:
        ExecutionStack.size()==2
        ExecutionStack.current() == s2
        ExecutionStack.withinWorkflow()

        when:
        ExecutionStack.push(s3)
        then:
        ExecutionStack.size()==3
        ExecutionStack.current() == s3
        ExecutionStack.withinWorkflow()

        when:
        def result = ExecutionStack.pop()
        then:
        result == s3
        ExecutionStack.size()==2
        ExecutionStack.current() == s2
        ExecutionStack.withinWorkflow()

        when:
        result = ExecutionStack.pop()
        then:
        result == s2
        ExecutionStack.size()==1
        ExecutionStack.current() == s1
        !ExecutionStack.withinWorkflow()

        when:
        result = ExecutionStack.pop()
        then:
        result == s1
        ExecutionStack.size()==0
        !ExecutionStack.withinWorkflow()

    }

}
