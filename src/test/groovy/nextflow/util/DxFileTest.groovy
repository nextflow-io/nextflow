package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DxFileTest extends Specification {

    def testEqualsAndHashCode() {

        when:
        def file1 = new DxFile('file-xxx', 'alpha')
        def file2 = new DxFile('file-xxx', 'alpha')
        def file3 = new DxFile('file-123', 'beta')

        then:
        file1 == file2
        file1 != file3

        file1.hashCode() == file2.hashCode()
        file2.hashCode() != file3.hashCode()

    }


    def testFileName() {

        when:
        def file = new DxFile('file-B7zqYX00j5831QBf1VV007xZ')

        then:
        file.name == 'split.nf'
        file.baseName == 'split'
        file.length() == 342
        file.size() == 342


    }

    def testText () {

        when:
        def file = new DxFile('file-B7zqYX00j5831QBf1VV007xZ')
        println file.text

        then:
        true

    }



}
