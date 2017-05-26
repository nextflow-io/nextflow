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
        def ch = Channel.from('A'..'Z')
        def sampler = new RandomSampleOp(ch, 10)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }


    def 'should produce random sample given a short channel' () {

        given:
        def ch = Channel.from('A'..'J')
        def sampler = new RandomSampleOp(ch, 20)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }

    def 'should produce random sample given a channel emitting the same number of items as the buffer' () {

        given:
        def ch = Channel.from('A'..'J')
        def sampler = new RandomSampleOp(ch, 10)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }

}
