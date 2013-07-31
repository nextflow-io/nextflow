package nextflow.executor

import com.dnanexus.DXJSON
import nextflow.exception.MissingFileException
import nextflow.util.DxFile
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DnaNexusExecutorTest extends Specification {

    def 'test outputs' () {

        setup:
        def outStr = '''
                   "output":{
                      "file1":"file-Z0C09V9F9WJQ7FJRKFKLLDFG",
                      "file2":"file-7F7G7G7S7HTKGLH0D9A8A8D8",
                      "file_x.fa":"file-B7vV85J5J6LLJ33K5LL6N4LY",
                      "file_y.fa":"file-B7vVF9G99GFMDK5KDNDJFY5Y",
                      "file_z.fa":"file-7FJFJ5J6SLA9FGKS6DHF8RFJ",
                      "exit_code":"0"
                   }
                '''

        def output = DXJSON.parseJson(outStr)


        // a single file is specification
        when:
        def file0 = DnaNexusExecutor.getFiles( output, 'file_z.fa' )
        then:
        file0 == new DxFile(name:'file_z.fa', id:'file-7FJFJ5J6SLA9FGKS6DHF8RFJ')


        // question mark wildcard is used
        when:
        def files1 = DnaNexusExecutor.getFiles( output, 'file?' )
        then:
        files1 == [
                new DxFile(name:'file_x.fa', id:'file-B7vV85J5J6LLJ33K5LL6N4LY'),
                new DxFile(name:'file_y.fa', id:'file-B7vVF9G99GFMDK5KDNDJFY5Y'),
                new DxFile(name:'file_z.fa', id:'file-7FJFJ5J6SLA9FGKS6DHF8RFJ')
        ]

        // star wildcard used, returns the matching list of DxFile(s)
        when:
        def files2 = DnaNexusExecutor.getFiles( output, 'file_*.fa' )
        then:
        files2 == [
                new DxFile(name:'file1', id:'file-Z0C09V9F9WJQ7FJRKFKLLDFG'),
                new DxFile(name:'file2', id:'file-7F7G7G7S7HTKGLH0D9A8A8D8')
        ]

        // missing file
        when:
        DnaNexusExecutor.getFiles(output,'missing.file')
        then:
        thrown(MissingFileException)

        // empty list
        when:
        def listEmpty = DnaNexusExecutor.getFiles(output, '*.gtf')
        then:
        listEmpty == []




    }

}
