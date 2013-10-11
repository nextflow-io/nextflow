package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BlankSeparatedListTest extends Specification {

    def testToString() {

        expect:
        new BlankSeparatedList('a'..'z').toString() == ('a'..'z').join(' ')
        "${new BlankSeparatedList('a'..'z')}" == ('a'..'z').join(' ')

    }
}
