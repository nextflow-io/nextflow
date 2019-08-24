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
