package nextflow.prov

import static test.TestHelper.*

import nextflow.config.ConfigParser
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import test.Dsl2Spec
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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
        def upstream = upstreamTasksOf('p2')
        upstream.size() == 1
        upstream.first.name == 'p1'
    }

    def 'should branch two process'() {

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
        def t1 = upstreamTasksOf('p2 (1)')
        t1.first.name == 'p1 (1)'
        t1.size() ==  1

        and:
        def t2 = upstreamTasksOf('p3 (1)')
        t2.first.name == 'p1 (2)'
        t2.size() ==  1

        and:
        def t3 = upstreamTasksOf('p3 (2)')
        t3.first.name == 'p1 (3)'
        t3.size() ==  1
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
        def t1 = upstreamTasksOf('p2 (1)')
        t1.first.name == 'p1 (1)'
        t1.size() ==  1
        
        and:
        def t2 = upstreamTasksOf('p2 (2)')
        t2.first.name == 'p1 (1)'
        t2.size() ==  1

        and:
        def t3 = upstreamTasksOf('p2 (3)')
        t3.first.name == 'p1 (2)'
        t3.size() ==  1

        and:
        def t4 = upstreamTasksOf('p2 (4)')
        t4.first.name == 'p1 (2)'
        t4.size() ==  1
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
        def t1 = upstreamTasksOf('p2')
        t1.name == ['p1 (1)', 'p1 (2)', 'p1 (3)']
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
        def t1 = upstreamTasksOf('p2 (1)')
        t1.name == ['p1 (1)', 'p1 (2)', 'p1 (3)']

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
        def t1 = upstreamTasksOf('p2')
        t1.name == ['p1 (1)', 'p1 (2)', 'p1 (3)']

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
        def t1 = upstreamTasksOf('p2')
        t1.name == ['p1 (1)', 'p1 (2)', 'p1 (3)']

    }

}
