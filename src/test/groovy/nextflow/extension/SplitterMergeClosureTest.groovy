package nextflow.extension

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitterMergeClosureTest extends Specification {

    def 'should merge splits' () {

        when:
        def merge = new SplitterMergeClosure([1,2])
        def output1 = ['pair_1', 'a', ['x','y'], 'any']
        def output2 = ['pair_1', ['a','b'], 'x', 'any']
        def result = merge.call( [output1, output2] as Object[] )
        then:
        result == ['pair_1', 'a', 'x', 'any']

        when:
        merge = new SplitterMergeClosure([2,3])
        output1 = ['pair_1', 'any' , 'a', ['x','y']]
        output2 = ['pair_1', 'any', ['a','b'], 'x']
        result = merge.call( [output1, output2] as Object[] )
        then:
        result == ['pair_1', 'any', 'a', 'x']

        when:
        merge = new SplitterMergeClosure([0,1])
        output1 = ['a', ['x','y'], 'pair_1', 'any' ]
        output2 = [['a','b'], 'x', 'pair_1', 'any']
        result = merge.call( [output1, output2] as Object[] )
        then:
        result == ['a', 'x', 'pair_1', 'any']
    }
}
