package nextflow.extension
import static test.TestHelper.gunzip

import java.nio.file.Files

import nextflow.Channel
import spock.lang.Specification
import spock.lang.Timeout
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class SplitFastqOperatorTest extends Specification {

    String READS = '''
        @SRR636272.19519409/1
        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
        +
        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
        @SRR636272.13995011/1
        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
        +
        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
        @SRR636272.21107783/1
        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
        +
        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
        @SRR636272.23331539/1
        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
        +
        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
        '''.stripIndent().leftTrim()

    String READS2 = '''
        @SRR636272.19519409/2
        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
        +
        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
        @SRR636272.13995011/2
        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
        +
        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
        @SRR636272.21107783/2
        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
        +
        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
        @SRR636272.23331539/2
        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
        +
        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
        '''.stripIndent().leftTrim()


    def 'should split a fastq' () {

        when:
        def target = Channel.from(READS).splitFastq(by:2)
        then:
        target.val == '''
            @SRR636272.19519409/1
            GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
            +
            CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
            @SRR636272.13995011/1
            GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
            +
            BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
            '''
            .stripIndent().leftTrim()

        target.val == '''
            @SRR636272.21107783/1
            CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
            +
            BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
            @SRR636272.23331539/1
            GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
            +
            CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
            '''
                .stripIndent().leftTrim()

        target.val == Channel.STOP
    }

    def 'should split a fastq to gzip chunks' () {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def target = Channel.from(READS).splitFastq(by:2, compress:true, file:folder)
        then:
        gunzip(target.val) == '''
            @SRR636272.19519409/1
            GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
            +
            CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
            @SRR636272.13995011/1
            GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
            +
            BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
            '''
                .stripIndent().leftTrim()

        gunzip(target.val) == '''
            @SRR636272.21107783/1
            CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
            +
            BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
            @SRR636272.23331539/1
            GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
            +
            CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
            '''
                .stripIndent().leftTrim()

        target.val == Channel.STOP

        cleanup:
        folder.deleteDir()
    }

    def 'should split read pairs' () {

        when:
        def result = Channel.from([['sample_id',READS,READS2]]).splitFastq(by:1, elem:[1,2]).toList().val

        then:
        result.size() ==4

        result[0][0] == 'sample_id'
        result[0][1] == '''
                        @SRR636272.19519409/1
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()
        result[0][2] == '''
                        @SRR636272.19519409/2
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()

        result[1][0] == 'sample_id'
        result[1][1] == '''
                        @SRR636272.13995011/1
                        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
                        +
                        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
                       '''.stripIndent().leftTrim()
        result[1][2] == '''
                        @SRR636272.13995011/2
                        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
                        +
                        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
                        '''.stripIndent().leftTrim()

        result[2][0] == 'sample_id'
        result[2][1] == '''
                        @SRR636272.21107783/1
                        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
                        +
                        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
                        '''.stripIndent().leftTrim()
        result[2][2] == '''
                        @SRR636272.21107783/2
                        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
                        +
                        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
                        '''.stripIndent().leftTrim()

        result[3][0] == 'sample_id'
        result[3][1] == '''
                        @SRR636272.23331539/1
                        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
                        +
                        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
                        '''.stripIndent().leftTrim()
        result[3][2] == '''
                        @SRR636272.23331539/2
                        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
                        +
                        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
                        '''.stripIndent().leftTrim()
    }

    def 'should split PE reads' () {
        given:
        def file1 = TestHelper.createInMemTempFile('one.fq', READS)
        def file2 = TestHelper.createInMemTempFile('two.fq', READS2)
        def channel
        def result

        when:
        channel = Channel.from([['sample_id',file1,file2]]).splitFastq(by:1, pe:true)
        result = channel.val
        then:
        result[0] == 'sample_id'
        result[1] == '''
                        @SRR636272.19519409/1
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()
        result[2] == '''
                        @SRR636272.19519409/2
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'sample_id'
        result[1] == '''
                        @SRR636272.13995011/1
                        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
                        +
                        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
                       '''.stripIndent().leftTrim()
        result[2] == '''
                        @SRR636272.13995011/2
                        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
                        +
                        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'sample_id'
        result[1] == '''
                        @SRR636272.21107783/1
                        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
                        +
                        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
                        '''.stripIndent().leftTrim()
        result[2] == '''
                        @SRR636272.21107783/2
                        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
                        +
                        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'sample_id'
        result[1] == '''
                        @SRR636272.23331539/1
                        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
                        +
                        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
                        '''.stripIndent().leftTrim()
        result[2] == '''
                        @SRR636272.23331539/2
                        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
                        +
                        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result == Channel.STOP

    }

    def 'split fastq to file chunks' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file_a_1 = TestHelper.createInMemTempFile('aaa_1.fq', READS)
        def file_a_2 = TestHelper.createInMemTempFile('aaa_2.fq', READS2)
        def file_b_1 = TestHelper.createInMemTempFile('bbb_1.fq', READS)
        def file_b_2 = TestHelper.createInMemTempFile('bbb_2.fq', READS2)
        def channel
        def result

        when:
        channel = Channel
                    .from([ ['aaa_id',file_a_1,file_a_2], ['bbb_id',file_b_1,file_b_2] ])
                    .splitFastq(by:1, pe:true, file:folder)
        result = channel.val
        then:
        result[0] == 'aaa_id'
        result[1].name == 'aaa_1.1.fq'
        result[2].name == 'aaa_2.1.fq'
        result[1].text == '''
                        @SRR636272.19519409/1
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()
        result[2].text == '''
                        @SRR636272.19519409/2
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'aaa_id'
        result[1].name == 'aaa_1.2.fq'
        result[2].name == 'aaa_2.2.fq'
        result[1].text == '''
                        @SRR636272.13995011/1
                        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
                        +
                        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
                       '''.stripIndent().leftTrim()
        result[2].text == '''
                        @SRR636272.13995011/2
                        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
                        +
                        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'aaa_id'
        result[1].name == 'aaa_1.3.fq'
        result[2].name == 'aaa_2.3.fq'
        result[1].text == '''
                        @SRR636272.21107783/1
                        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
                        +
                        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
                        '''.stripIndent().leftTrim()
        result[2].text == '''
                        @SRR636272.21107783/2
                        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
                        +
                        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'aaa_id'
        result[1].name == 'aaa_1.4.fq'
        result[2].name == 'aaa_2.4.fq'
        result[1].text == '''
                        @SRR636272.23331539/1
                        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
                        +
                        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
                        '''.stripIndent().leftTrim()
        result[2].text == '''
                        @SRR636272.23331539/2
                        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
                        +
                        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'bbb_id'
        result[1].name == 'bbb_1.1.fq'
        result[2].name == 'bbb_2.1.fq'
        result[1].text == '''
                        @SRR636272.19519409/1
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()
        result[2].text == '''
                        @SRR636272.19519409/2
                        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
                        +
                        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
                        '''.stripIndent().leftTrim()

        when:
        result = channel.val
        then:
        result[0] == 'bbb_id'
        result[1].name == 'bbb_1.2.fq'
        result[2].name == 'bbb_2.2.fq'

        when:
        result = channel.val
        then:
        result[0] == 'bbb_id'
        result[1].name == 'bbb_1.3.fq'
        result[2].name == 'bbb_2.3.fq'

        when:
        result = channel.val
        then:
        result[0] == 'bbb_id'
        result[1].name == 'bbb_1.4.fq'
        result[2].name == 'bbb_2.4.fq'

        when:
        result = channel.val
        then:
        result == Channel.STOP
        // finally check the existence of the
        folder.resolve('.chunks.aaa_1.fq').exists()
        folder.resolve('.chunks.aaa_2.fq').exists()
        folder.resolve('.chunks.bbb_1.fq').exists()
        folder.resolve('.chunks.bbb_2.fq').exists()
    }
}
