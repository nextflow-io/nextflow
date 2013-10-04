package nextflow.executor
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.node.ObjectNode
import groovy.io.FileType
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DnaNexusExecutorTest extends Specification {




    protected ObjectNode makeDXLink(String objectId) {
        return DXJSON.getObjectBuilder().put('$dnanexus_link', objectId).build();
    }



    def void testFindFiles() {

        setup:
        def folder = Paths.get(URI.create('dxfs:///test_1/'))
        Files.createDirectories( folder )

        def target1 = folder.resolve('test_file1.txt')
        def target2 = folder.resolve('test_file2.txt')
        def target3 = folder.resolve('diff_name.fa')
        Files.copy(new ByteArrayInputStream("Hello1".getBytes()), target1);
        Files.copy(new ByteArrayInputStream("Hello2".getBytes()), target2);
        Files.copy(new ByteArrayInputStream("Hello3".getBytes()), target3);

        when:
        Path folderUnderTest = Paths.get(URI.create('dxfs:///test_1/'))

        List files = []
        folderUnderTest.eachFileMatch(FileType.FILES, ~/test_.*/ ) { Path it -> files << it}

        then:
        files.size() == 2


        cleanup:
        folderUnderTest?.deleteDir()

    }


}