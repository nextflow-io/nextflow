/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.processor

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Global
import nextflow.ISession
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.file.FileHolder
import nextflow.file.FilePorter
import nextflow.script.ScriptType
import nextflow.util.ArrayBag
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskInputResolverTest extends Specification {


    def 'should return single item or collection'() {

        setup:
        def path1 = Paths.get('file1')
        def path2 = Paths.get('file2')
        def path3 = Paths.get('file3')

        when:
        def list = [ FileHolder.get(path1, 'x_file_1') ]
        def result = TaskInputResolver.singleItemOrList(list, true, ScriptType.SCRIPTLET)
        then:
        result.toString() == 'x_file_1'

        when:
        list = [ FileHolder.get(path1, 'x_file_1') ]
        result = TaskInputResolver.singleItemOrList(list, false, ScriptType.SCRIPTLET)
        then:
        result*.toString() == ['x_file_1']

        when:
        list = [ FileHolder.get(path1, 'x_file_1'), FileHolder.get(path2, 'x_file_2'), FileHolder.get(path3, 'x_file_3') ]
        result = TaskInputResolver.singleItemOrList(list, false, ScriptType.SCRIPTLET)
        then:
        result*.toString() == [ 'x_file_1',  'x_file_2',  'x_file_3']

    }


    def 'should expand wildcards'() {

        /*
         * The name do not contain any wildcards *BUT* when multiple files are provide
         * an index number is added to the specified name
         */
        when:
        def list1 = TaskInputResolver.expandWildcards('file_name', [FileHolder.get('x')])
        def list2 = TaskInputResolver.expandWildcards('file_name', [FileHolder.get('x'), FileHolder.get('y')] )
        then:
        list1 *. stageName  == ['file_name']
        list2 *. stageName  == ['file_name1', 'file_name2']


        /*
         * The star wildcard: when a single item is provided, it is simply ignored
         * When a collection of files is provided, the name is expanded to the index number
         */
        when:
        list1 = TaskInputResolver.expandWildcards('file*.fa', [FileHolder.get('x')])
        list2 = TaskInputResolver.expandWildcards('file_*.fa', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['file.fa']
        list2 *. stageName == ['file_1.fa', 'file_2.fa', 'file_3.fa']

        /*
         * The question mark wildcards *always* expand to an index number
         */
        when:
        def p0 = [FileHolder.get('0')]
        def p1_p4 = (1..4).collect { FileHolder.get(it.toString()) }
        def p1_p12 = (1..12).collect { FileHolder.get(it.toString()) }
        list1 = TaskInputResolver.expandWildcards('file?.fa', p0 )
        list2 = TaskInputResolver.expandWildcards('file_???.fa', p1_p4 )
        def list3 = TaskInputResolver.expandWildcards('file_?.fa', p1_p12 )
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list3 instanceof ArrayBag
        list1 *. stageName == ['file1.fa']
        list2 *. stageName == ['file_001.fa', 'file_002.fa', 'file_003.fa', 'file_004.fa']
        list3 *. stageName == ['file_1.fa', 'file_2.fa', 'file_3.fa', 'file_4.fa', 'file_5.fa', 'file_6.fa', 'file_7.fa', 'file_8.fa', 'file_9.fa', 'file_10.fa', 'file_11.fa', 'file_12.fa']

        when:
        list1 = TaskInputResolver.expandWildcards('*', [FileHolder.get('a')])
        list2 = TaskInputResolver.expandWildcards('*', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['a']
        list2 *. stageName == ['x','y','z']

        when:
        list1 = TaskInputResolver.expandWildcards('dir1/*', [FileHolder.get('a')])
        list2 = TaskInputResolver.expandWildcards('dir2/*', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['dir1/a']
        list2 *. stageName == ['dir2/x','dir2/y','dir2/z']

        when:
        list1 = TaskInputResolver.expandWildcards('/dir/file*.fa', [FileHolder.get('x')])
        list2 = TaskInputResolver.expandWildcards('dir/file_*.fa', [FileHolder.get('x'), FileHolder.get('y'), FileHolder.get('z')])
        then:
        list1 instanceof ArrayBag
        list2 instanceof ArrayBag
        list1 *. stageName == ['dir/file.fa']
        list2 *. stageName == ['dir/file_1.fa', 'dir/file_2.fa', 'dir/file_3.fa']

        when:
        list1 = TaskInputResolver.expandWildcards('dir/*', [FileHolder.get('file.fa')])
        list2 = TaskInputResolver.expandWildcards('dir/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/file.fa']
        list2 *. stageName == ['dir/titi.fa', 'dir/toto.fa']

        when:
        list1 = TaskInputResolver.expandWildcards('dir/*/*', [FileHolder.get('file.fa')])
        list2 = TaskInputResolver.expandWildcards('dir/*/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/1/file.fa']
        list2 *. stageName == ['dir/1/titi.fa', 'dir/2/toto.fa']

        when:
        list1 = TaskInputResolver.expandWildcards('dir/foo*/*', [FileHolder.get('file.fa')])
        list2 = TaskInputResolver.expandWildcards('dir/foo*/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/foo1/file.fa']
        list2 *. stageName == ['dir/foo1/titi.fa', 'dir/foo2/toto.fa']

        when:
        list1 = TaskInputResolver.expandWildcards('dir/??/*', [FileHolder.get('file.fa')])
        list2 = TaskInputResolver.expandWildcards('dir/??/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/01/file.fa']
        list2 *. stageName == ['dir/01/titi.fa', 'dir/02/toto.fa']

        when:
        list1 = TaskInputResolver.expandWildcards('dir/bar??/*', [FileHolder.get('file.fa')])
        list2 = TaskInputResolver.expandWildcards('dir/bar??/*', [FileHolder.get('titi.fa'), FileHolder.get('file.fq', 'toto.fa')])
        then:
        list1 *. stageName == ['dir/bar01/file.fa']
        list2 *. stageName == ['dir/bar01/titi.fa', 'dir/bar02/toto.fa']
    }

    @Unroll
    def 'should expand wildcards rule' () {

        expect:
        TaskInputResolver.expandWildcards0(pattern, 'stage-name.txt', index, size ) == expected

        where:
        pattern             | index | size  | expected
        // just wildcard
        '*'                 | 1     | 1     | 'stage-name.txt'
        '*'                 | 1     | 10    | 'stage-name.txt'
        // wildcard on the file name and single item in the collection
        'foo.txt'           | 1     | 1     | 'foo.txt'
        'foo*.fa'           | 1     | 1     | 'foo.fa'
        'foo?.fa'           | 1     | 1     | 'foo1.fa'
        'foo??.fa'          | 1     | 1     | 'foo01.fa'
        // wildcard on the file name and many items in the collection
        'foo*.fa'           | 1     | 3     | 'foo1.fa'
        'foo*.fa'           | 3     | 3     | 'foo3.fa'
        'foo?.fa'           | 1     | 3     | 'foo1.fa'
        'foo?.fa'           | 3     | 3     | 'foo3.fa'
        'foo??.fa'          | 1     | 3     | 'foo01.fa'
        'foo??.fa'          | 3     | 3     | 'foo03.fa'
        // wildcard on parent path
        'dir/*/foo.txt'     | 1     | 1     | 'dir/1/foo.txt'
        'dir/foo*/bar.txt'  | 1     | 1     | 'dir/foo1/bar.txt'
        'dir/foo?/bar.txt'  | 2     | 2     | 'dir/foo2/bar.txt'
        'dir/foo??/bar.txt' | 2     | 2     | 'dir/foo02/bar.txt'
        // wildcard on parent path and name
        'dir/*/'            | 1     | 1     | 'dir/1/stage-name.txt'
        'dir/*/*'           | 1     | 1     | 'dir/1/stage-name.txt'
        'dir/*/*'           | 1     | 10    | 'dir/1/stage-name.txt'
        'dir/*/foo*.txt'    | 1     | 1     | 'dir/1/foo.txt'
        'dir/*/foo*.txt'    | 1     | 2     | 'dir/1/foo1.txt'
        'dir/*/foo?.txt'    | 2     | 2     | 'dir/2/foo2.txt'
        'dir/???/foo?.txt'  | 5     | 10    | 'dir/005/foo5.txt'
    }

    @Unroll
    def 'should replace question marks' () {
        expect:
        TaskInputResolver.replaceQuestionMarkWildcards(pattern, index) == expected

        where:
        pattern         | index | expected
        'foo.txt'       | 1     | 'foo.txt'
        'foo?.txt'      | 1     | 'foo1.txt'
        'foo???.txt'    | 2     | 'foo002.txt'
        'foo?_???.txt'  | 3     | 'foo3_003.txt'
        'foo??.txt'     | 9999  | 'foo9999.txt'

    }

    def "should return a file holder" () {

        given:
        FileHolder holder
        def tempFolder = Files.createTempDirectory('test')
        def localFile = Files.createTempFile(tempFolder, 'test','test')
        Global.session = Mock(ISession)
        Global.session.workDir >> tempFolder

        /*
         * when the input file is on the local file system
         * simple return a reference to it in the holder object
         */
        when:
        holder = TaskInputResolver.normalizeInputToFile(localFile,null)
        then:
        holder.sourceObj == localFile
        holder.storePath == localFile.toRealPath()
        holder.stageName == localFile.getFileName().toString()

        /*
         * any generic input that is not a file is converted to a string
         * and save to the local file system
         */
        when:
        holder = TaskInputResolver.normalizeInputToFile("text data string",'simple_file_name.txt')
        then:
        holder.sourceObj == "text data string"
        holder.storePath.fileSystem == FileSystems.default
        holder.storePath.text == "text data string"
        holder.stageName == 'simple_file_name.txt'

        cleanup:
        tempFolder?.deleteDir()
    }

    def 'should normalise to path' () {
        expect:
        TaskInputResolver.normalizeToPath('/foo/bar') == '/foo/bar' as Path
        and:
        TaskInputResolver.normalizeToPath('file:///foo/bar') == '/foo/bar' as Path
        and:
        TaskInputResolver.normalizeToPath(Paths.get('foo.txt')) == Paths.get('foo.txt')

        when:
        TaskInputResolver.normalizeToPath('abc')
        then:
        thrown(ProcessUnrecoverableException)

        when:
        TaskInputResolver.normalizeToPath(null)
        then:
        thrown(ProcessUnrecoverableException)
    }

    def 'should normalize files' () {
        given:
        def batch = Mock(FilePorter.Batch)
        def executor = Mock(Executor)
        def PATH = Paths.get('/some/path')
        def resolver = new TaskInputResolver(Mock(TaskRun), batch, executor)

        when:
        def result = resolver.normalizeInputToFiles(PATH.toString(), 0, true)
        then:
        1 * executor.isForeignFile(PATH) >> false
        0 * batch.addToForeign(PATH) >> null
        and:
        result.size() == 1
        result[0] == new FileHolder(PATH)

    }

}
