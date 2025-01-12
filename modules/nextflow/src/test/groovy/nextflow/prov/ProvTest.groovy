package nextflow.prov

import static test.TestHelper.*

import nextflow.config.ConfigParser
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import spock.lang.Ignore
import spock.lang.Timeout
import test.Dsl2Spec
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ProvTest extends Dsl2Spec {

    def setup() {
        Prov.clear()
        TaskId.clear()
        TaskProcessor.allTasks.clear()
    }

    ConfigObject globalConfig() {
        new ConfigParser().parse('''
        process.fair = true
        ''')
    }

    def 'should chain two process'() {

        when:
        dsl_eval(globalConfig(), '''
            workflow {
                p1 | map { x-> x } | map { x-> x+1 }  | p2 
            }
            
            process p1 { 
              output: val(x) 
              exec: 
                x =1 
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1']
    }

    def 'should track provenance with branch operator'() {

        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,10,20) \
                    | p1 \
                    | branch { left: it <=10; right: it >10 } \
                    | set { result }
                
                result.left | p2
                result.right | p3                                
            }
            
            process p1 {
              input: val(x) 
              output: val(y) 
              exec: 
                y = x+1 
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
            
            process p3 {
              input: val(x)
              exec: 
                println x
            }
        ''')
        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)']

        and:
         upstreamTasksOf('p3 (1)')
                .name == ['p1 (2)']

        and:
        upstreamTasksOf('p3 (2)')
                .name == ['p1 (3)']
    }

    def 'should track provenance with flatMap operator' () {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2) \
                    | p1 \
                    | flatMap \
                    | p2
            }
            
            process p1 {
              input: val(x) 
              output: val(y) 
              exec: 
                y = [x, x*x] 
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')
        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)']

        and:
        upstreamTasksOf('p2 (2)')
                .name == ['p1 (1)']

        and:
        upstreamTasksOf('p2 (3)')
                .name == ['p1 (2)']

        and:
        upstreamTasksOf('p2 (4)')
                .name == ['p1 (2)']
    }

    def 'should track the provenance of two processes and reduce operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3) \
                    | p1 \
                    | reduce {a,b -> return a+b} \
                    | p2
            }
            
            process p1 {
              input: val(x) 
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (1)', 'p1 (2)', 'p1 (3)']
    }

    def 'should track the provenance of two tasks and collectFile operator' () {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of('a','b','c') \
                    | p1 \
                    | collectFile(name: 'sample.txt') \
                    | p2
            }
            
            process p1 {
              input: val(x) 
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: file(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)', 'p1 (2)', 'p1 (3)']

    }

    def 'should track the provenance of two tasks and toList operator' () {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of('a','b','c') | p1 | toList | p2
            }
            
            process p1 {
              input: val(x) 
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (1)', 'p1 (2)', 'p1 (3)']

    }

    def 'should track the provenance of two tasks and toSortedList operator' () {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of('a','b','c') | p1 | toList | p2
            }
            
            process p1 {
              input: val(x) 
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (1)', 'p1 (2)', 'p1 (3)']

    }

    def 'should track provenance with filter operator'() {

        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,3,2) | p1 | filter { x-> x<30 } | p2
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x*10 
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)']
        then:
        upstreamTasksOf('p2 (2)')
                .name == ['p1 (3)']
    }

    def 'should track provenance with unique operator'() {

        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,2,3) | p1 | unique | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)']

        then:
        upstreamTasksOf('p2 (2)')
                .name == ['p1 (2)']

        then:
        upstreamTasksOf('p2 (3)')
                .name == ['p1 (4)']
    }


    def 'should track provenance with distinc operator'() {

        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,2,2,3) | p1 | unique | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        def upstream1 = upstreamTasksOf('p2 (1)')
        upstream1.size() == 1
        upstream1.first.name == 'p1 (1)'

        then:
        upstreamTasksOf('p2 (2)')
                .name == ['p1 (2)']

        then:
        upstreamTasksOf('p2 (3)')
                .name == ['p1 (5)']
    }

    def 'should track provenance with first operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3) | p1 | first | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (1)']
    }

    def 'should track provenance with take operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3,4,5) | p1 | take(2) | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)']
        then:
        upstreamTasksOf('p2 (2)')
                .name == ['p1 (2)']

    }

    def 'should track provenance with last operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3,4,5) | p1 | last | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (5)']
    }

    def 'should track provenance with collect operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3) | p1 | collect | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (1)', 'p1 (2)', 'p1 (3)']
    }

    def 'should track provenance with count value operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.value(1) | p1 | count | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1']
    }

    def 'should track provenance with count many operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3) | p1 | count | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (1)', 'p1 (2)', 'p1 (3)']
    }

    @Ignore // this should be review
    def 'should track provenance with min operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(3,2,1) | p1 | min | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2')
                .name == ['p1 (3)']
    }

    def 'should track provenance with buffer operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                channel.of(1,2,3,4,5) | p1 | buffer(size:2, remainder:true) | p2 
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p2 (1)')
                .name == ['p1 (1)', 'p1 (2)']
        and:
        upstreamTasksOf('p2 (2)')
                .name == ['p1 (3)', 'p1 (4)']

        and:
        upstreamTasksOf('p2 (3)')
                .name == ['p1 (5)']
    }

    @Ignore
    def 'should track provenance with mix operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                def c1 = channel.of(1,2) | p1 
                def c2 = channel.of(3,4) | p2 
                p1.out | mix(p2.out) | p3  
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p3 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p3 (1)')
                .name == ['p1 (1)']

    }

    def 'should track provenance with join operator'() {
        when:
        dsl_eval(globalConfig(), '''
            workflow {
                def c1 = channel.of(['a',10], ['b',20]) | p1 
                def c2 = channel.of(['a',11], ['b',21], ['c',31]) | p2 
                p1.out | join(p2.out, remainder:true) | p3
            }
            
            process p1 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p2 { 
              input: val(x)
              output: val(y) 
              exec: 
                y = x
            }
            
            process p3 {
              input: val(x)
              exec: 
                println x
            }
        ''')

        then:
        upstreamTasksOf('p3 (1)')
                .name.sort() == ['p1 (1)', 'p2 (1)']

        and:
        upstreamTasksOf('p3 (2)')
            .name.sort() == ['p1 (2)', 'p2 (2)']

        and:
        upstreamTasksOf('p3 (3)')
            .name.sort() == ['p2 (3)']

    }
}
