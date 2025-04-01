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
 *
 */

package nextflow.data.cid.fs

import nextflow.data.cid.model.WorkflowResults
import nextflow.data.cid.serde.CidEncoder
import nextflow.file.FileHelper
import nextflow.serde.gson.GsonEncoder

import java.nio.file.Files
import java.time.Instant

import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll

/**
 * CID Path Tests
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CidPathTest extends Specification {

    @Shared def wdir = Files.createTempDirectory('wdir')
    @Shared def cid = wdir.resolve('.meta')
    @Shared def data = wdir.resolve('work')
    @Shared def fs = Mock(CidFileSystem)

    def cleanupSpec(){
        wdir.deleteDir()
    }

    def 'should create from URI' () {
        when:
        def path = new CidPath(fs, new URI( URI_STRING ))
        then:
        path.filePath == PATH
        path.fragment == FRAGMENT
        path.query == QUERY

        where:
        URI_STRING                      | PATH              | QUERY         | FRAGMENT
        "cid://1234/hola"               | "1234/hola"       | null          | null
        "cid://1234/hola#frag.sub"      | "1234/hola"       | null          | "frag.sub"
        "cid://1234/#frag.sub"          | "1234"            | null          | "frag.sub"
        "cid://1234/?q=a&b=c"           | "1234"            | "q=a&b=c"     | null
        "cid://1234/?q=a&b=c#frag.sub"  | "1234"            | "q=a&b=c"     | "frag.sub"
        "cid:///"                       | "/"               | null          | null
    }

    def 'should create correct cid Path' () {
        when:
            def cid = new CidPath(FS, PATH, MORE)
        then:
            cid.filePath == EXPECTED_FILE
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

        def cidFs = new FileHelper().getOrCreateFileSystemFor('cid', [enabled: true, store: [location: cid.parent.toString()]] )

        cid.resolve('12345/output1').mkdirs()
        cid.resolve('12345/path/to/file2.txt').mkdirs()
        cid.resolve('12345/.data.json').text = '{"type":"TaskRun"}'
        cid.resolve('12345/output1/.data.json').text = '{"type":"TaskOutput", "path": "' + outputFolder.toString() + '"}'
        cid.resolve('12345/path/to/file2.txt/.data.json').text = '{"type":"TaskOutput", "path": "' + outputFile.toString() + '"}'
        def time = Instant.now().toString()
        def wfResultsMetadata = new CidEncoder().withPrettyPrint(true).encode(new WorkflowResults(time, "cid://1234", [a: "cid://1234/a.txt"]))
        cid.resolve('5678/').mkdirs()
        cid.resolve('5678/.data.json').text = wfResultsMetadata

        expect: 'Get real path when CidPath is the output data or a subfolder'
        new CidPath(cidFs,'12345/output1' ).getTargetPath() == outputFolder
        new CidPath(cidFs,'12345/output1/some/path' ).getTargetPath() == outputSubFolder
        new CidPath(cidFs,'12345/output1/some/path/file1.txt').getTargetPath().text == outputSubFolderFile.text
        new CidPath(cidFs, '12345/path/to/file2.txt').getTargetPath().text == outputFile.text

        when: 'CidPath fs is null'
        new CidPath(null, '12345').getTargetPath()
        then:
        thrown(IllegalArgumentException)

        when: 'CidPath is empty'
        new CidPath(cidFs, '/').getTargetPath()
        then:
        thrown(IllegalArgumentException)

        when: 'CidPath is not an output data description'
        new CidPath(cidFs, '12345').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'CidPath is not subfolder of an output data description'
        new CidPath(cidFs, '12345/path').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'Cid does not exist'
        new CidPath(cidFs, '23456').getTargetPath()
        then:
        thrown(FileNotFoundException)

        when: 'Cid description'
        def result = new CidPath(cidFs, '5678').getTargetPath(true)
        then:
        result instanceof CidResultsPath
        result.text == wfResultsMetadata

        when: 'Cid description subobject'
        def result2 = new CidPath(cidFs, '5678#outputs').getTargetPath(true)
        then:
        result2 instanceof CidResultsPath
        result2.text == new GsonEncoder<Object>(){}.withPrettyPrint(true).encode([a: "cid://1234/a.txt"])

        when: 'Cid subobject does not exist'
        new CidPath(cidFs, '23456#notexists').getTargetPath(true)
        then:
        thrown(FileNotFoundException)

        cleanup:
        cid.resolve('12345').deleteDir()

    }

    def 'should get file name' () {
        when:
        def cid1 = new CidPath(fs, '1234567890/this/file.bam')
        then:
        cid1.getFileName() == new CidPath(null, 'file.bam')
    }

    def 'should get file parent' () {
        when:
        def cid1 = new CidPath(fs, '1234567890/this/file.bam')
        then:
        cid1.getParent() == new CidPath(fs, '1234567890/this')
        cid1.getParent().getParent() == new CidPath(fs, '1234567890')
        cid1.getParent().getParent().getParent() == new CidPath(fs, "/")
        cid1.getParent().getParent().getParent().getParent() == null
    }

    @Unroll
    def 'should get name count' () {
        expect:
        new CidPath(fs, PATH).getNameCount() == EXPECTED
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
        new CidPath(fs, PATH).getName(INDEX) == EXPECTED
        where:
        PATH        | INDEX | EXPECTED
        '123'       | 0     | new CidPath(fs, '123')
        '123/a'     | 1     | new CidPath(null, 'a')
        '123/a/'    | 1     | new CidPath(null, 'a')
        '123/a/b'   | 2     | new CidPath(null, 'b')
    }

    @Unroll
    def 'should get subpath' () {
        expect:
        new CidPath(fs, PATH).subpath(BEGIN,END) == EXPECTED
        where:
        PATH        | BEGIN | END   | EXPECTED
        '123'       | 0     | 1     | new CidPath(fs, '123')
        '123/a'     | 0     | 2     | new CidPath(fs, '123/a')
        '123/a/'    | 0     | 2     | new CidPath(fs, '123/a')
        '123/a'     | 1     | 2     | new CidPath(null, 'a')
        '123/a/'    | 1     | 2     | new CidPath(null, 'a')
        '123/a/b'   | 2     | 3     | new CidPath(null, 'b')
        '123/a/b'   | 1     | 3     | new CidPath(null, 'a/b')
    }

    def 'should normalize a path' () {
        expect:
        new CidPath(fs, '123').normalize() == new CidPath(fs, '123')
        new CidPath(fs, '123/a/b').normalize() == new CidPath(fs, '123/a/b')
        new CidPath(fs, '123/./a/b').normalize() == new CidPath(fs, '123/a/b')
        new CidPath(fs, '123/a/../a/b').normalize() == new CidPath(fs, '123/a/b')
    }

    @Unroll
    def 'should validate startWith' () {
        expect:
        new CidPath(fs,PATH).startsWith(OTHER) == EXPECTED
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
        new CidPath(fs,PATH).endsWith(OTHER) == EXPECTED
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
        new CidPath(fs,'1234/a/b/c').isAbsolute()
        new CidPath(fs,'1234/a/b/c').getRoot().isAbsolute()
        new CidPath(fs,'1234/a/b/c').getParent().isAbsolute()
        new CidPath(fs,'1234/a/b/c').normalize().isAbsolute()
        new CidPath(fs,'1234/a/b/c').getName(0).isAbsolute()
        new CidPath(fs,'1234/a/b/c').subpath(0,2).isAbsolute()
        and:
        !new CidPath(fs,'1234/a/b/c').getFileName().isAbsolute()
        !new CidPath(fs,'1234/a/b/c').getName(1).isAbsolute()
        !new CidPath(fs,'1234/a/b/c').subpath(1,3).isAbsolute()
    }

    @Unroll
    def 'should get root path' () {
        expect:
        new CidPath(fs,PATH).getRoot() == new CidPath(fs,EXPECTED)
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
        new CidPath(fs, '/')            | new CidPath(fs, '123/a/b/c')      | new CidPath(null, '123/a/b/c')
        new CidPath(fs,'123/a/')        | new CidPath(fs, '123/a/b/c')      | new CidPath(null, 'b/c')
        new CidPath(fs,'123/a/')        | new CidPath(fs, '321/a/')         | new CidPath(null, '../../321/a')
        new CidPath(null,'123/a')       | new CidPath(null, '123/a/b/c')    | new CidPath(null, 'b/c')
        new CidPath(null,'123/a')       | new CidPath(null, '321/a')        | new CidPath(null, '../../321/a')
        new CidPath(fs,'../a/')         | new CidPath(fs, '321/a')          | new CidPath(null, '../321/a')
        new CidPath(fs,'321/a/')        | new CidPath(fs, '../a')           | new CidPath(null, '../../a')
        new CidPath(null,'321/a/')      | new CidPath(null, '../a')         | new CidPath(null, '../../../a')
    }

    def 'relativize should throw exception' () {
        given:
        def cid1 = new CidPath(fs,'123/a/')
        def cid2 = new CidPath(null,'123/a/')
        def cid3 = new CidPath(null, '../a/b')
        when: 'comparing relative with absolute'
            cid1.relativize(cid2)
        then:
            thrown(IllegalArgumentException)

        when: 'undefined base path'
            cid3.relativize(cid2)
        then:
            thrown(IllegalArgumentException)
    }

    def 'should resolve path' () {
        when:
        def cid1 = new CidPath(fs, '123/a/b/c')
        def cid2 = new CidPath(fs, '321/x/y/z')
        def rel1 = new CidPath(null, 'foo')
        def rel2 = new CidPath(null, 'bar/')
        
        then:
        cid1.resolve(cid2) == cid2
        cid2.resolve(cid1) == cid1
        and:
        cid1.resolve(rel1) == new CidPath(fs,'123/a/b/c/foo')
        cid1.resolve(rel2) == new CidPath(fs,'123/a/b/c/bar')
        and:
        rel1.resolve(rel2) == new CidPath(null, 'foo/bar')
        rel2.resolve(rel1) == new CidPath(null, 'bar/foo')
    }

    def 'should resolve path as string' () {
        given:
        def pr = Mock(CidFileSystemProvider)
        def cidfs = Mock(CidFileSystem){
            provider() >> pr}


        def cid1 = new CidPath(cidfs, '123/a/b/c')

        expect:
        cid1.resolve('x/y') == new CidPath(cidfs, '123/a/b/c/x/y')
        cid1.resolve('/x/y/') == new CidPath(cidfs, '123/a/b/c/x/y')

        when:
        def result = cid1.resolve('cid://321')
        then:
        pr.getPath(CidPath.asUri('cid://321')) >> new CidPath(cidfs, '321')
        and:
        result == new CidPath(cidfs, '321')
    }
    
    @Unroll
    def 'should get to uri string' () {
        expect:
        new CidPath(fs, PATH).toUriString() == EXPECTED
        where:
        PATH            | EXPECTED
        '/'             | 'cid:///'
        '1234'         | 'cid://1234'
        '1234/a/b/c'   | 'cid://1234/a/b/c'
        ''             | 'cid:///'
    }

    @Unroll
    def 'should get string' () {
        expect:
        new CidPath(fs, PATH).toString() == EXPECTED
        where:
        PATH            | EXPECTED
        '/'             | '/'
        '1234'          | '1234'
        '1234/a/b/c'    | '1234/a/b/c'
        ''              | '/'
    }

    def 'should validate asString method'() {
        expect:
        CidPath.asUriString(FIRST, MORE as String[]) == EXPECTED

        where:
        FIRST       | MORE          | EXPECTED
        'foo'       | []            | 'cid://foo'
        'foo/'      | []            | 'cid://foo'
        '/foo'      | []            | 'cid://foo'
        and:
        'a'       | ['/b/']         | 'cid://a/b'
        'a'       | ['/b','c']      | 'cid://a/b/c'
        'a'       | ['/b','//c']    | 'cid://a/b/c'
        'a'       | ['/b/c', 'd']   | 'cid://a/b/c/d'
        '/a/'     | ['/b/c', 'd']   | 'cid://a/b/c/d'
    }
}
