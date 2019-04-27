package nextflow.script

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelArrayListTest extends Specification {

    def 'should get values' () {

        when:
        def arr = new ChannelArrayList()
        then:
        arr.fifth == null
        arr.second == null

        when:
        arr = new ChannelArrayList([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        then:
        arr.first == 1
        arr.second == 2
        arr.third == 3
        arr.fourth == 4
        arr.fifth == 5
        arr.sixth == 6
        arr.seventh == 7
        arr.eighth == 8
        arr.ninth == 9
        arr.tenth == 10

    }

    def 'should validate output spread' () {
        given:
        def out1 = new ChannelArrayList( ['a', 'b'] )
        def out2 = new ChannelArrayList( ['x', 'y', 'z'])

        expect:
        ChannelArrayList.spread([1, 2, 3]) \
            == [1, 2, 3]

        ChannelArrayList.spread([out1]) \
            == ['a', 'b']

        ChannelArrayList.spread([out1, 'p', 'q', out2] ) \
            == ['a', 'b', 'p', 'q', 'x', 'y', 'z']
    }
}
