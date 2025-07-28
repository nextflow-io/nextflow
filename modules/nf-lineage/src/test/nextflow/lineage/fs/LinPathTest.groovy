/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.fs

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.time.OffsetDateTime

import nextflow.lineage.LinUtils
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowOutput
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.util.CacheHelper
import org.junit.Rule
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll
import test.OutputCapture
/**
 * LID Path Tests
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinPathTest extends Specification {

    @Shared
    Path wdir
    @Shared
    Path data
    @Shared
    def fs = Mock(LinFileSystem)

    def setupSpec(){
        wdir = Files.createTempDirectory("wdir")
        data = wdir.resolve('work')
    }

    def cleanupSpec(){
        wdir.deleteDir()
    }

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should create from URI' () {
        when:
        def path = new LinPath(fs, new URI( URI_STRING ))
        then:
        path.filePath == PATH
        path.fragment == FRAGMENT
        path.query == QUERY

        where:
        URI_STRING                                  | PATH              | QUERY         | FRAGMENT
        "lid://1234/hola"                           | "1234/hola"       | null          | null
        "lid://1234/hola#workflow.repository"       | "1234/hola"       | null          | "workflow.repository"
        "lid://1234/#workflow.repository"           | "1234"            | null          | "workflow.repository"
        "lid://1234/?q=a&b=c"                       | "1234"            | "q=a&b=c"     | null
        "lid://1234/?q=a&b=c#workflow.repository"   | "1234"            | "q=a&b=c"     | "workflow.repository"
        "lid:///"                                   | "/"               | null          | null
    }

    def 'should throw exception if fragment contains an unknown property'() {
        when:
        new LinPath(fs, new URI ("lid://1234/hola#no-exist"))
        then:
        thrown(IllegalArgumentException)

    }

    def 'should warn if query is specified'() {
        when:
        new LinPath(fs, new URI("lid://1234/hola?query"))
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0].endsWith("Query string is not supported for Lineage URI: `lid://1234/hola?query` -- it will be ignored")
    }

    def 'should create correct lid Path' () {
        when:
            def lid = new LinPath(FS, PATH, MORE)
        then:
            lid.filePath == EXPECTED_FILE
        where:
        FS      | PATH                  | MORE                      | EXPECTED_FILE
        fs      | '/'                   | [] as String[]            | '/'
        null    | '/'                   | [] as String[]            | '/'
        fs      | '/'                   | ['a','b'] as String[]     | 'a/b'
        null    | '/'                   | ['a','b'] as String[]     | 'a/b'
        fs      | ''                    | [] as String[]            | '/'
        null    | ''                    | [] as String[]            | '/'
        fs      | ''                    | ['a','b'] as String[]     | 'a/b'
        null    | ''                    | ['a','b'] as String[]     | 'a/b'
        fs      | '1234'                | [] as String[]            | '1234'
        null    | '1234'                | [] as String[]            | '1234'
        fs      | '1234'                | ['a','b'] as String[]     | '1234/a/b'
        null    | '1234'                | ['a','b'] as String[]     | '1234/a/b'
        fs      | '1234/c'              | [] as String[]            | '1234/c'
        null    | '1234/c'              | [] as String[]            | '1234/c'
        fs      | '1234/c'              | ['a','b'] as String[]     | '1234/c/a/b'
        null    | '1234/c'              | ['a','b'] as String[]     | '1234/c/a/b'
        fs      | '/1234/c'             | [] as String[]            | '1234/c'
        null    | '/1234/c'             | [] as String[]            | '1234/c'
        fs      | '/1234/c'             | ['a','b'] as String[]     | '1234/c/a/b'
        null    | '/1234/c'             | ['a','b'] as String[]     | '1234/c/a/b'
        fs      | '../c'                | ['a','b'] as String[]     | 'c/a/b'
        null    | '../c'                | ['a','b'] as String[]     | '../c/a/b'
        fs      | '../c'                | [] as String[]            | 'c'
        null    | '../c'                | [] as String[]            | '../c'
        fs      | '..'                  | [] as String[]            | '/'
        null    | '..'                  | [] as String[]            | '..'
        fs      | '/..'                 | [] as String[]            | '/'
        null    | '/..'                 | [] as String[]            | '/'
        fs      | './1234/c'            | ['a','b'] as String[]     | '1234/c/a/b'
        null    | './1234/c'            | ['a','b'] as String[]     | '1234/c/a/b'
        fs      | './1234/c'            | [] as String[]            | '1234/c'
        null    | './1234/c'            | [] as String[]            | '1234/c'
        fs      | '1234'                | ['/'] as String[]         | '1234'
        null    | '1234'                | ['/'] as String[]         | '1234'
        null    | '../../a/b'           | [] as String[]            | '../../a/b'
        fs      | '1234/'               | [] as String[]            | '1234'
        null    | '1234/'               | [] as String[]            | '1234'
    }

    def 'should get target path' () {
        given:
        def outputFolder = data.resolve('output')
        def outputSubFolder = outputFolder.resolve('some/path')
        outputSubFolder.mkdirs()
        def outputSubFolderFile = outputSubFolder.resolve('file1.txt')
        outputSubFolderFile.text = "this is file1"
        def outputFile = data.resolve('file2.txt')
        outputFile.text = "this is file2"

        def lidFs = new LinFileSystemProvider().newFileSystem(new URI("lid:///"), [enabled: true, store: [location: wdir.toString()]])

        wdir.resolve('12345/output1').mkdirs()
        wdir.resolve('12345/path/to/file2.txt').mkdirs()
        wdir.resolve('12345/.data.json').text = '{"version":"lineage/v1beta1","kind":"TaskRun","spec":{"name":"test"}}'
        wdir.resolve('12345/output1/.data.json').text = '{"version":"lineage/v1beta1","kind":"FileOutput","spec":{"path": "' + outputFolder.toString() + '"}}'
        wdir.resolve('12345/path/to/file2.txt/.data.json').text = '{"version":"lineage/v1beta1","kind":"FileOutput","spec":{"path": "' + outputFile.toString() + '"}}'
        def time = OffsetDateTime.now()
        def wfResultsMetadata = new LinEncoder().withPrettyPrint(true).encode(new WorkflowOutput(time, "lid://1234", [new Parameter( "Path", "a", "lid://1234/a.txt")]))
        wdir.resolve('5678/').mkdirs()
        wdir.resolve('5678/.data.json').text = wfResultsMetadata

        expect: 'Get real path when LinPath is the output data or a subfolder'
        new LinPath(lidFs, '12345/output1').getTargetPath() == outputFolder
        new LinPath(lidFs,'12345/output1/some/path').getTargetPath() == outputSubFolder
        new LinPath(lidFs,'12345/output1/some/path/file1.txt').getTargetPath().text == outputSubFolderFile.text
        new LinPath(lidFs, '12345/path/to/file2.txt').getTargetPath().text == outputFile.text

        when: 'LinPath fs is null'
        new LinPath(null, '12345').getTargetPath()
        then:
        thrown(IllegalArgumentException)

        when: 'LinPath is empty'
        new LinPath(lidFs, '/').getTargetPath()
        then:
        thrown(IllegalArgumentException)

        when: 'LinPath is not an output data description'
        new LinPath(lidFs, '12345').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'LinPath is not subfolder of an output data description'
        new LinPath(lidFs, '12345/path').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'LinPath subfolder of an output data description does not exist'
        new LinPath(lidFs, '12345/output1/other/path').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'Lid does not exist'
        new LinPath(lidFs, '23456').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'LinPath is not an output data description but is a workflow/task run'
        def result = new LinPath(lidFs, '12345').getTargetOrIntermediatePath()
        then:
        result instanceof LinIntermediatePath

        when: 'LinPath is a subpath of a workflow/task run data description'
        result = new LinPath(lidFs, '12345/path').getTargetOrIntermediatePath()
        then:
        result instanceof LinIntermediatePath

        when: 'Lid description'
        result = new LinPath(lidFs, '5678').getTargetOrMetadataPath()
        then:
        result instanceof LinMetadataPath
        result.text == wfResultsMetadata

        when: 'Lid description subobject'
        def result2 = new LinPath(lidFs, '5678#output').getTargetOrMetadataPath()
        then:
        result2 instanceof LinMetadataPath
        result2.text == LinUtils.encodeSearchOutputs([new Parameter("Path","a", "lid://1234/a.txt")], true)

        when: 'Lid subobject does not exist'
        new LinPath(lidFs, '23456#notexists').getTargetOrMetadataPath()
        then:
        thrown(IllegalArgumentException)
    }

    def 'should get subobjects as path' (){
        given:
        def lidFs = new LinFileSystemProvider().newFileSystem(new URI("lid:///"), [enabled: true, store: [location: wdir.toString()]])
        def wf = new WorkflowRun(new Workflow([],"repo", "commit"), "sessionId", "runId", [new Parameter("String", "param1", "value1")])

        when: 'workflow repo in workflow run'
        Path p = LinPath.getMetadataAsTargetPath(wf, lidFs, "123456", "workflow.repository")
        then:
        p instanceof LinMetadataPath
        p.text == '"repo"'

        when: 'outputs'
        def outputs = new WorkflowOutput(OffsetDateTime.now(), "lid://123456", [new Parameter("Collection", "samples", ["sample1", "sample2"])])
        lidFs.store.save("123456#output", outputs)
        Path p2 = LinPath.getMetadataAsTargetPath(wf, lidFs, "123456", "output")
        then:
        p2 instanceof LinMetadataPath
        p2.text == LinUtils.encodeSearchOutputs([new Parameter("Collection", "samples", ["sample1", "sample2"])], true)

        when: 'child does not exists'
        LinPath.getMetadataAsTargetPath(wf, lidFs, "123456", "no-exist")
        then:
        def exception = thrown(FileNotFoundException)
        exception.message == "Target path '123456#no-exist' does not exist"

        when: 'outputs does not exists'
        LinPath.getMetadataAsTargetPath(wf, lidFs, "6789", "output")
        then:
        def exception1 = thrown(FileNotFoundException)
        exception1.message == "Target path '6789#output' does not exist"

        when: 'null object'
        LinPath.getMetadataAsTargetPath(null, lidFs, "123456", "no-exist")
        then:
        def exception2 = thrown(FileNotFoundException)
        exception2.message == "Target path '123456' does not exist"

        cleanup:
        wdir.resolve("123456").deleteDir()
    }

    def 'should get file name' () {
        expect:
        new LinPath(fs, PATH).getFileName() == EXPECTED
        where:
        PATH                        | EXPECTED
        '1234567890/this/file.bam'  | new LinPath(null, 'file.bam')
        '12345/hola?query#output'   | new LinPath("query", "output", "hola", null)

    }

    def 'should get file parent' () {
        when:
        def lid1 = new LinPath(fs, '1234567890/this/file.bam')
        then:
        lid1.getParent() == new LinPath(fs, '1234567890/this')
        lid1.getParent().getParent() == new LinPath(fs, '1234567890')
        lid1.getParent().getParent().getParent() == new LinPath(fs, "/")
        lid1.getParent().getParent().getParent().getParent() == null
    }

    @Unroll
    def 'should get name count' () {
        expect:
        new LinPath(fs, PATH).getNameCount() == EXPECTED
        where:
        PATH        | EXPECTED
        '/'         | 0
        '123'       | 1
        '123/a'     | 2
        '123/a/'    | 2
        '123/a/b'   | 3
        ''          | 0
    }

    @Unroll
    def 'should get name by index' () {
        expect:
        new LinPath(fs, PATH).getName(INDEX) == EXPECTED
        where:
        PATH                | INDEX | EXPECTED
        '123'               | 0     | new LinPath(fs, '123')
        '123/a'             | 1     | new LinPath(null, 'a')
        '123/a/'            | 1     | new LinPath(null, 'a')
        '123/a/b'           | 2     | new LinPath(null, 'b')
        '123/a?q#output'    | 1     | new LinPath(null, 'a?q#output')
    }

    @Unroll
    def 'should get subpath' () {
        expect:
        new LinPath(fs, PATH).subpath(BEGIN,END) == EXPECTED
        where:
        PATH        | BEGIN | END   | EXPECTED
        '123'       | 0     | 1     | new LinPath(fs, '123')
        '123/a'     | 0     | 2     | new LinPath(fs, '123/a')
        '123/a/'    | 0     | 2     | new LinPath(fs, '123/a')
        '123/a'     | 1     | 2     | new LinPath(null, 'a')
        '123/a/'    | 1     | 2     | new LinPath(null, 'a')
        '123/a/b'   | 2     | 3     | new LinPath(null, 'b')
        '123/a/b'   | 1     | 3     | new LinPath(null, 'a/b')
    }

    def 'should normalize a path' () {
        expect:
        new LinPath(fs, '123').normalize() == new LinPath(fs, '123')
        new LinPath(fs, '123/a/b').normalize() == new LinPath(fs, '123/a/b')
        new LinPath(fs, '123/./a/b').normalize() == new LinPath(fs, '123/a/b')
        new LinPath(fs, '123/a/../a/b').normalize() == new LinPath(fs, '123/a/b')
    }

    @Unroll
    def 'should validate startWith' () {
        expect:
        new LinPath(fs,PATH).startsWith(OTHER) == EXPECTED
        where:
        PATH            | OTHER         | EXPECTED
        '12345/a/b'     | '12345'       | true
        '12345/a/b'     | '12345/a'     | true
        '12345/a/b'     | '12345/a/b'   | true
        and:
        '12345/a/b'     | '12345/b'     | false
        '12345/a/b'     | 'xyz'         | false
    }

    @Unroll
    def 'should validate endsWith' () {
        expect:
        new LinPath(fs,PATH).endsWith(OTHER) == EXPECTED
        where:
        PATH            | OTHER         | EXPECTED
        '12345/a/b'     | 'b'           | true
        '12345/a/b'     | 'a/b'         | true
        '12345/a/b'     | '12345/a/b'   | true
        and:
        '12345/a/b'     | '12345/b'     | false
        '12345/a/b'     | 'xyz'         | false
    }

    def 'should validate isAbsolute' () {
        expect:
        new LinPath(fs,'1234/a/b/c').isAbsolute()
        new LinPath(fs,'1234/a/b/c').getRoot().isAbsolute()
        new LinPath(fs,'1234/a/b/c').getParent().isAbsolute()
        new LinPath(fs,'1234/a/b/c').normalize().isAbsolute()
        new LinPath(fs,'1234/a/b/c').getName(0).isAbsolute()
        new LinPath(fs,'1234/a/b/c').subpath(0,2).isAbsolute()
        and:
        !new LinPath(fs,'1234/a/b/c').getFileName().isAbsolute()
        !new LinPath(fs,'1234/a/b/c').getName(1).isAbsolute()
        !new LinPath(fs,'1234/a/b/c').subpath(1,3).isAbsolute()
    }

    @Unroll
    def 'should get root path' () {
        expect:
        new LinPath(fs,PATH).getRoot() == new LinPath(fs,EXPECTED)
        where:
        PATH            | EXPECTED
        '12345'         | '/'
        '12345/a'       | '/'
    }

    def 'should relativize path' () {
        expect:
        BASE_PATH.relativize(PATH) == EXPECTED
        where :
        BASE_PATH                       | PATH                              | EXPECTED
        new LinPath(fs, '/')       | new LinPath(fs, '123/a/b/c')   | new LinPath(null, '123/a/b/c')
        new LinPath(fs,'123/a/')   | new LinPath(fs, '123/a/b/c')   | new LinPath(null, 'b/c')
        new LinPath(fs,'123/a/')   | new LinPath(fs, '321/a/')      | new LinPath(null, '../../321/a')
        new LinPath(null,'123/a')  | new LinPath(null, '123/a/b/c') | new LinPath(null, 'b/c')
        new LinPath(null,'123/a')  | new LinPath(null, '321/a')     | new LinPath(null, '../../321/a')
        new LinPath(fs,'../a/')    | new LinPath(fs, '321/a')       | new LinPath(null, '../321/a')
        new LinPath(fs,'321/a/')   | new LinPath(fs, '../a')        | new LinPath(null, '../../a')
        new LinPath(null,'321/a/') | new LinPath(null, '../a')      | new LinPath(null, '../../../a')
    }

    def 'relativize should throw exception' () {
        given:
        def lid1 = new LinPath(fs,'123/a/')
        def lid2 = new LinPath(null,'123/a/')
        def lid3 = new LinPath(null, '../a/b')
        when: 'comparing relative with absolute'
            lid1.relativize(lid2)
        then:
            thrown(IllegalArgumentException)

        when: 'undefined base path'
            lid3.relativize(lid2)
        then:
            thrown(IllegalArgumentException)
    }

    def 'should resolve path' () {
        when:
        def lid1 = new LinPath(fs, '123/a/b/c')
        def lid2 = new LinPath(fs, '321/x/y/z')
        def rel1 = new LinPath(null, 'foo')
        def rel2 = new LinPath(null, 'bar/')
        
        then:
        lid1.resolve(lid2) == lid2
        lid2.resolve(lid1) == lid1
        and:
        lid1.resolve(rel1) == new LinPath(fs,'123/a/b/c/foo')
        lid1.resolve(rel2) == new LinPath(fs,'123/a/b/c/bar')
        and:
        rel1.resolve(rel2) == new LinPath(null, 'foo/bar')
        rel2.resolve(rel1) == new LinPath(null, 'bar/foo')
    }

    def 'should resolve path as string' () {
        given:
        def pr = Mock(LinFileSystemProvider)
        def lidfs = Mock(LinFileSystem){
            provider() >> pr}


        def lid1 = new LinPath(lidfs, '123/a/b/c')

        expect:
        lid1.resolve('x/y') == new LinPath(lidfs, '123/a/b/c/x/y')
        lid1.resolve('/x/y/') == new LinPath(lidfs, '123/a/b/c/x/y')

        when:
        def result = lid1.resolve('lid://321')
        then:
        pr.getPath(LinPath.asUri('lid://321')) >> new LinPath(lidfs, '321')
        and:
        result == new LinPath(lidfs, '321')
    }

    def 'should throw illegal exception when not correct scheme' (){
        when: 'creation'
        new LinPath(fs, new URI("http://1234"))
        then:
        thrown(IllegalArgumentException)

        when: 'asUri'
        LinPath.asUri("http://1234")
        then:
        thrown(IllegalArgumentException)

        when: 'asUri'
        LinPath.asUri("")
        then:
        thrown(IllegalArgumentException)

    }

    def 'should throw provider mismatch exception when different path types' () {
        given:
        def pr = Mock(LinFileSystemProvider)
        def fs = Mock(LinFileSystem){
            provider() >> pr}
        and:
        def lid = new LinPath(fs, '123/a/b/c')

        when: 'resolve with path'
        lid.resolve(Path.of('d'))
        then:
        thrown(ProviderMismatchException)

        when: 'resolve with uri string'
        lid.resolve(Path.of('http://1234'))
        then:
        thrown(ProviderMismatchException)

        when: 'relativize'
        lid.relativize(Path.of('d'))
        then:
        thrown(ProviderMismatchException)
    }

    def 'should throw exception for unsupported methods' () {
        given:
        def pr = Mock(LinFileSystemProvider)
        def fs = Mock(LinFileSystem){
            provider() >> pr}
        and:
        def lid = new LinPath(fs, '123/a/b/c')

        when: 'to file'
        lid.toFile()
        then:
        thrown(UnsupportedOperationException)

        when: 'register'
        lid.register(null, null,null)
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should throw exception for incorrect index'() {
        when: 'getting name with negative index'
        new LinPath(fs, "1234").getName(-1)
        then:
        thrown(IllegalArgumentException)

        when: 'getting name with larger index tha namecount'
        new LinPath(fs, "1234").getName(2)
        then:
        thrown(IllegalArgumentException)

        when: 'getting subpath with negative index'
        new LinPath(fs, "1234").subpath(-1,1)
        then:
        thrown(IllegalArgumentException)

        when: 'getting subpath with larger index tha namecount'
        new LinPath(fs, "1234").subpath(0,2)
        then:
        thrown(IllegalArgumentException)

    }

    @Unroll
    def 'should get to uri string' () {
        expect:
        new LinPath(fs, PATH).toUriString() == EXPECTED
        where:
        PATH            | EXPECTED
        '/'             | 'lid:///'
        '1234'         | 'lid://1234'
        '1234/a/b/c'   | 'lid://1234/a/b/c'
        ''             | 'lid:///'
    }

    @Unroll
    def 'should get string' () {
        expect:
        new LinPath(fs, PATH).toString() == EXPECTED
        where:
        PATH            | EXPECTED
        '/'             | '/'
        '1234'          | '1234'
        '1234/a/b/c'    | '1234/a/b/c'
        ''              | '/'
    }

    @Unroll
    def 'should validate asString method'() {
        expect:
        LinPath.asUriString(FIRST, MORE as String[]) == EXPECTED

        where:
        FIRST       | MORE          | EXPECTED
        'foo'       | []            | 'lid://foo'
        'foo/'      | []            | 'lid://foo'
        '/foo'      | []            | 'lid://foo'
        and:
        'a'       | ['/b/']         | 'lid://a/b'
        'a'       | ['/b','c']      | 'lid://a/b/c'
        'a'       | ['/b','//c']    | 'lid://a/b/c'
        'a'       | ['/b/c', 'd']   | 'lid://a/b/c/d'
        '/a/'     | ['/b/c', 'd']   | 'lid://a/b/c/d'
    }

    @Unroll
    def 'should check is lid uri string' () {
        expect:
        LinPath.isLidUri(STR) == EXPECTED

        where:
        STR             | EXPECTED
        null            | false
        ''              | false
        'foo'           | false
        '/foo'          | false
        'lid:/foo'      | false
        'lid:foo'       | false
        'lid/foo'       | false
        and:
        'lid://'        | true
        'lid:///'       | true
        'lid://foo/bar' | true
    }

    def 'should detect equals'(){
        expect:
        new LinPath(FS1, PATH1).equals(new LinPath(FS2, PATH2)) == EXPECTED
        where:
        FS1            | FS2            | PATH1             | PATH2             | EXPECTED
        null           | fs             | "12345/path"      | "12345/path"      | false
        fs             | null           | "12345/path"      | "12345/path"      | false
        null           | null           | "12345/"          | "12345/path"      | false
        fs             | fs             | "12345/"          | "12345/path"      | false
        and:
        null           | null           | "12345/path"      | "12345/path"      | true
        fs             | fs             | "12345/path"      | "12345/path"      | true
        null           | null           | "12345/"          | "12345"           | true
        fs             | fs             | "12345/"          | "12345 "          | true
    }

    def 'should validate correct hash'(){
        given:
        def file = wdir.resolve("file.txt")
        file.text = "this is a data file"
        def hash = CacheHelper.hasher(file).hash().toString()
        def correctData = new FileOutput(file.toString(), new Checksum(hash,"nextflow", "standard"))
        when:
        def error = LinPath.validateDataOutput(correctData)
        then:
        !error

        cleanup:
        file.delete()
    }

    def 'should warn with incorrect hash'(){
        given:
        def file = wdir.resolve("file.txt")
        file.text = "this is a data file"
        def hash = CacheHelper.hasher(file).hash().toString()
        def correctData = new FileOutput(file.toString(), new Checksum("abscd","nextflow", "standard"))
        when:
        def error = LinPath.validateDataOutput(correctData)
        then:
        error == "Checksum of '$file' does not match with lineage metadata"

        cleanup:
        file.delete()
    }

    def 'should warn when hash algorithm is not supported'(){
        given:
        def file = wdir.resolve("file.txt")
        file.text = "this is a data file"
        def hash = CacheHelper.hasher(file).hash().toString()
        def correctData = new FileOutput(file.toString(), new Checksum(hash,"not-supported", "standard"))
        when:
        def error = LinPath.validateDataOutput(correctData)
        then:
        error == "Checksum of '$file' can't be validated - algorithm 'not-supported' is not supported"

        cleanup:
        file.delete()
    }

    def 'should throw exception when file not found validating hash'(){
        when:
        def correctData = new FileOutput("not/existing/file", new Checksum("120741","nextflow", "standard"))
        LinPath.validateDataOutput(correctData)

        then:
        thrown(FileNotFoundException)
    }

    def 'should validate path' () {
        given:
        def outputFolder = data.resolve('output')
        outputFolder.mkdirs()
        def outputFile = outputFolder.resolve('file1.txt')
        outputFile.text = "this is file1"
       def encoder = new LinEncoder()
        def hash = CacheHelper.hasher(outputFile).hash().toString()
        def correctData = new FileOutput(outputFile.toString(), new Checksum(hash,"nextflow", "standard"))
        def incorrectData = new FileOutput(outputFile.toString(), new Checksum("incorrectHash","nextflow", "standard"))
        wdir.resolve('12345/output/file1.txt').mkdirs()
        wdir.resolve('12345/output/file2.txt').mkdirs()
        wdir.resolve('12345/output/file1.txt/.data.json').text = encoder.encode(correctData)
        wdir.resolve('12345/output/file2.txt/.data.json').text = encoder.encode(incorrectData)
        def lidFs = new LinFileSystemProvider().newFileSystem(new URI("lid:///"), [enabled: true, store: [location: wdir.toString()]])

        when:
        def succeed = new LinPath(lidFs, '12345/output/file1.txt').validate()
        then:
        succeed

        when:
        succeed = new LinPath(lidFs, '12345/output/file2.txt').validate()
        then:
        !succeed

        when:
        succeed = new LinPath(lidFs, '12345/output/file3.txt').validate()
        then:
        !succeed

        cleanup:
        outputFile.delete()
    }


}
