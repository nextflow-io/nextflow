package nextflow.executor

import nextflow.script.FileInParam
import nextflow.processor.FileHolder
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractExecutorTest extends Specification {




    def testStagingFilesScript() {
        setup:
        def owner = Mock(Script)
        def executor = [:] as AbstractExecutor

        def param1 = new FileInParam(owner, 'file.txt')
        def param2 = new FileInParam(owner, 'seq_*.fa')
        Map<FileInParam,List<FileHolder>> files = [:]
        files[param1] = [FileHolder.get('/home/data/sequences', 'file.txt')]
        files[param2] = [FileHolder.get('/home/data/file1','seq_1.fa'), FileHolder.get('/home/data/file2','seq_2.fa'), FileHolder.get('/home/data/file3','seq_3.fa') ]

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


        when:
        files = [:]
        files[param1] = [FileHolder.get('/data/file', 'file.txt')]
        files[param2] = [FileHolder.get('/data/seq','seq_1.fa') ]
        script = executor.stagingFilesScript(files, '; ')
        lines = script.readLines()

        then:
        lines[0] == 'rm -f file.txt; rm -f seq_1.fa; ln -s /data/file file.txt; ln -s /data/seq seq_1.fa; '
        lines.size() == 1

    }


}
