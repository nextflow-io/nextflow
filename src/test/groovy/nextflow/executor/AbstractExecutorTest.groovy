package nextflow.executor
import java.nio.file.Paths

import nextflow.processor.FileInParam
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractExecutorTest extends Specification {



    def testWildcardExpand() {

        setup:
        def executor = [:] as AbstractExecutor

        /*
         * The name do not contain any wildcards *BUT* when multiple files are provide
         * an index number is added to the specified name
         */
        when:
        def list1 = executor.expandWildcards('file_name', Paths.get('x'))
        def list2 = executor.expandWildcards('file_name', [Paths.get('x'), Paths.get('y')])
        then:
        list1 == ['file_name']
        list2 == ['file_name1', 'file_name2']


        /*
         * The star wildcard: when a single item is provided, it is simply ignored
         * When a collection of files is provided, the name is expanded to the index number
         */
        when:
        list1 = executor.expandWildcards('file*.fa', Paths.get('x'))
        list2 = executor.expandWildcards('file_*.fa', [Paths.get('x'), Paths.get('y'), Paths.get('z')])
        then:
        list1 == ['file.fa']
        list2 == ['file_1.fa', 'file_2.fa', 'file_3.fa']

        /*
         * The question mark wildcards *always* expand to an index number
         */
        when:
        list1 = executor.expandWildcards('file?.fa', Paths.get('0'))
        list2 = executor.expandWildcards('file_???.fa', 1..4 )
        def list3 = executor.expandWildcards('file_?.fa', 1..12 )
        then:
        list1 == ['file1.fa']
        list2 == ['file_001.fa', 'file_002.fa', 'file_003.fa', 'file_004.fa']
        list3 == ['file_1.fa', 'file_2.fa', 'file_3.fa', 'file_4.fa', 'file_5.fa', 'file_6.fa', 'file_7.fa', 'file_8.fa', 'file_9.fa', 'file_10.fa', 'file_11.fa', 'file_12.fa']

        when:
        list1 = executor.expandWildcards('*', Paths.get('a'))
        list2 = executor.expandWildcards('*', [Paths.get('x'), Paths.get('y'), Paths.get('z')])
        then:
        list1 == ['a']
        list2 == ['x','y','z']

        when:
        executor.expandWildcards('*', 0)
        then:
        thrown(IllegalArgumentException)

    }


    def testStagingFilesScript() {
        setup:
        def owner = Mock(Script)
        def executor = [:] as AbstractExecutor

        def param1 = new FileInParam(owner, 'file.txt')
        def param2 = new FileInParam(owner, 'seq_*.fa')
        Map<FileInParam,Object> files = [:]
        files[param1] = Paths.get('/home/data/sequences')
        files[param2] = [Paths.get('/home/data/file1'), Paths.get('/home/data/file2'), Paths.get('/home/data/file3') ]

        when:
        def script = executor.stagingFilesScript(files)
        def lines = script.readLines()

        then:
        lines[0] == 'rm -f file.txt'
        lines[1] == 'rm -f seq_1.fa'
        lines[2] == 'rm -f seq_2.fa'
        lines[3] == 'rm -f seq_3.fa'
        lines[4] == 'ln -s /home/data/sequences file.txt'
        lines[5] == 'ln -s /home/data/file1 seq_1.fa'
        lines[6] == 'ln -s /home/data/file2 seq_2.fa'
        lines[7] == 'ln -s /home/data/file3 seq_3.fa'
        lines.size() == 8

    }


}
