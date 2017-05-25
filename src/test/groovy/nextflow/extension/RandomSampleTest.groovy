package nextflow.extension

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RandomSampleTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should produce random sample' () {

        given:
        def ch = Channel.from(0,1,2,3,4,5,6,7,8,9)
        def sampler = new RandomSampleOp(ch, 5)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 5
        result.unique().size() == 5 
        result != [0,1,2,3,4]
        result[0] in 0..9
        result[1] in 0..9
        result[2] in 0..9
        result[3] in 0..9
        result[4] in 0..9
    }


    def 'should produce random sample given a short channel' () {

        given:
        def ch = Channel.from(0,1,2,3,4)
        def sampler = new RandomSampleOp(ch, 10)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 5
        result.unique().size() == 5
        result != [0,1,2,3,4]
        result[0] in 0..9
        result[1] in 0..9
        result[2] in 0..9
        result[3] in 0..9
        result[4] in 0..9
    }

    def 'should produce random sample given a channel emitting the same number of items as the buffer' () {

        given:
        def ch = Channel.from(0,1,2,3,4)
        def sampler = new RandomSampleOp(ch, 5)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 5
        result.unique().size() == 5
        result != [0,1,2,3,4]
        result[0] in 0..9
        result[1] in 0..9
        result[2] in 0..9
        result[3] in 0..9
        result[4] in 0..9
    }

}
