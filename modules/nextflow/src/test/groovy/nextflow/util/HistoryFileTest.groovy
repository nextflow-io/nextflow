/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.util
import static nextflow.util.HistoryFile.Record

import java.nio.file.Files

import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HistoryFileTest extends Specification {

    static String FILE_TEXT = '''
b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98\tnextflow run examples/ampa.nf --in data/sample.fa
b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98\tnextflow run examples/ampa.nf --in data/sample.fa -resume
58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8\tnextflow run examples/basic.nf --in data/sample.fa
2016-07-24 16:43:16\t-\tevil_pike\tOK\t6b9515aba6\te710da1b-ce06-482f-bbcf-987a507f85d1\t.nextflow run hello
2016-07-24 16:43:34\t-\tgigantic_keller\tOK\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello
2016-07-24 16:43:34\t-\tsmall_cirum\tOK\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello -resume
2016-07-25 09:58:01\t5 min\tmodest_bartik\tERR\t6b9515aba6\t5910a50f-8656-4765-aa79-f07cef912062\t.nextflow run hello
'''

    def 'test add and get and find' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()

        when:
        def history = new HistoryFile(file)
        then:
        history.getLast() == null

        when:
        def id1 = UUID.randomUUID()
        def id2 = UUID.randomUUID()
        def id3 = UUID.randomUUID()
        def now = System.currentTimeMillis()
        def d1 = new Date(now - 50_000)
        def d2 = new Date(now - 30_000)
        def d3 = new Date(now - 10_000)
        history.write( 'hello_world', id1, 'abc', [1,2,3], d1 )
        history.write( 'super_star', id2, '123', [1,2,3], d2 )
        history.write( 'slow_food', id3, 'xyz', [1,2,3], d3 )

        then:
        history.getLast() == new Record(sessionId: id3, runName: 'slow_food', timestamp: d3, command: '1 2 3')
        history.checkExistsById( id1.toString() )
        history.checkExistsById( id2.toString() )
        history.checkExistsById( id3.toString() )
        !history.checkExistsById( UUID.randomUUID().toString() )

    }


    def 'should return a session ID given a short version of it' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT

        when:
        def history = new HistoryFile(file)
        then:
        history.findById('b8a3c4cf') == [new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98')]
        history.findById('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98') ==  [new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98')]
        history.findById('58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8') == [new Record('58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8')]
        history.findById('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9') == [new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9','gigantic_keller'), new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9','small_cirum')]
        history.findById('5910a50f') == [new Record('5910a50f-8656-4765-aa79-f07cef912062','modest_bartik')]
        history.findById('5910a50x') == []

        history.checkExistsById('5910a50f')
        !history.checkExistsById('5910a50x')

        when:
        history.findById('5')
        then:
        def e = thrown(AbortOperationException)
        e.message == '''
                Which session ID do you mean?
                    58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8
                    5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9
                    5910a50f-8656-4765-aa79-f07cef912062
                '''
                .stripIndent().leftTrim()

    }


    def 'should return a session ID given a run name' () {

        given:
        def file = Files.createTempFile('test',null)
        file.text = FILE_TEXT

        when:
        def history = new HistoryFile(file)
        then:
        history.getByName('lazy_pike') == null
        history.getByName('evil_pike') == new Record('e710da1b-ce06-482f-bbcf-987a507f85d1', 'evil_pike')
        history.getByName('gigantic_keller') == new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9', 'gigantic_keller')

        cleanup:
        file?.delete()
    }

    def 'should verify uuid char' () {
        expect:
        HistoryFile.isUuidChar('-' as char)
        HistoryFile.isUuidChar('0' as char)
        HistoryFile.isUuidChar('3' as char)
        HistoryFile.isUuidChar('9' as char)
        HistoryFile.isUuidChar('a' as char)
        HistoryFile.isUuidChar('b' as char)
        HistoryFile.isUuidChar('f' as char)
        !HistoryFile.isUuidChar('q' as char)
        !HistoryFile.isUuidChar('!' as char)
    }

    def 'should verify uuid string' () {

        expect:
        HistoryFile.isUuidString('b')
        HistoryFile.isUuidString('b8a3c4cf')
        HistoryFile.isUuidString('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98')

        !HistoryFile.isUuidString('hello_world')
    }

    def 'should find entries before the specified one' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.findBefore('5a6d3877') == [
                new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98'),
                new Record('58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8'),
                new Record('e710da1b-ce06-482f-bbcf-987a507f85d1', 'evil_pike')
        ]

        history.findBefore('unknown') == []

    }

    def 'should find entries after the specified one' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.findAfter('evil_pike') == [
                new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9', 'gigantic_keller'),
                new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9', 'small_cirum'),
                new Record('5910a50f-8656-4765-aa79-f07cef912062', 'modest_bartik'),
        ]


        history.findAfter('unknown') == []
    }


    def 'should return all IDs' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.findAll() == [
                new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98'),
                new Record('58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8'),
                new Record('e710da1b-ce06-482f-bbcf-987a507f85d1', 'evil_pike'),
                new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9', 'gigantic_keller'),
                new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9', 'small_cirum'),
                new Record('5910a50f-8656-4765-aa79-f07cef912062', 'modest_bartik')
        ]
    }

    def 'should find by a session by a name of ID' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.findByIdOrName('last') == [new Record('5910a50f-8656-4765-aa79-f07cef912062', 'modest_bartik')]
        history.findByIdOrName('evil_pike') == [new Record('e710da1b-ce06-482f-bbcf-987a507f85d1', 'evil_pike')]
        history.findByIdOrName('b8a3c4cf') ==  [new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98')]
        history.findByIdOrName('5a6d3877') == [new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9','gigantic_keller'), new Record('5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9','small_cirum')]
        history.findByIdOrName('unknown') == []
    }

    def 'should delete history entry ' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        when:
        history.deleteEntry( new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98') )
        history.deleteEntry( new Record('e710da1b-ce06-482f-bbcf-987a507f85d1', 'evil_pike') )
        then:
        history.text == '''
                58d8dd16-ce77-4507-ba1a-ec1ccc9bd2e8\tnextflow run examples/basic.nf --in data/sample.fa
                2016-07-24 16:43:34\t-\tgigantic_keller\tOK\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello
                2016-07-24 16:43:34\t-\tsmall_cirum\tOK\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello -resume
                2016-07-25 09:58:01\t5 min\tmodest_bartik\tERR\t6b9515aba6\t5910a50f-8656-4765-aa79-f07cef912062\t.nextflow run hello
                '''
                .stripIndent()
    }

    def 'should list entries by id ' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.findById('b8a3c4cf') == [ new Record('b8a3c4cf-17e4-49c6-a4cf-4fd8ddbeef98') ]
        history.findById('e710da1b') == [ new Record('e710da1b-ce06-482f-bbcf-987a507f85d1', 'evil_pike') ]
        history.findById('e710dadd') == []
    }


    def 'should return true when find an entry by name' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.checkExistsByName('evil_pike')
        !history.checkExistsByName('missing')
    }


    def 'should return unique name' () {

        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        when:
        def name = history.generateNextName()
        then:
        name != null
        !history.findAllRunNames().contains(name)

    }

    def 'should return all run names in the history file' () {
        given:
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = FILE_TEXT
        def history = new HistoryFile(file)

        expect:
        history.findAllRunNames() == ['evil_pike', 'gigantic_keller', 'small_cirum', 'modest_bartik'] as Set
    }

    def 'should update the history entries ' () {

        given:
        def source = '''
2016-07-24 16:43:16\t-\tevil_pike\t-\t6b9515aba6\te710da1b-ce06-482f-bbcf-987a507f85d1\t.nextflow run hello
2016-07-24 16:43:34\t-\tgigantic_keller\t-\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello
2016-07-24 16:43:34\t-\tsmall_cirum\t-\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello -resume
2016-07-25 09:58:01\t5 min\tmodest_bartik\tERR\t6b9515aba6\t5910a50f-8656-4765-aa79-f07cef912062\t.nextflow run hello
'''
        def file = Files.createTempFile('test',null)
        file.deleteOnExit()
        file.text = source
        def history = new HistoryFile(file)


        when:
        def when = HistoryFile.TIMESTAMP_FMT.parse('2016-07-24 16:53:16')
        history.update('evil_pike',true,when)
        then:
        history.text == '''
2016-07-24 16:43:16\t10m\tevil_pike\tOK\t6b9515aba6\te710da1b-ce06-482f-bbcf-987a507f85d1\t.nextflow run hello
2016-07-24 16:43:34\t-\tgigantic_keller\t-\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello
2016-07-24 16:43:34\t-\tsmall_cirum\t-\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello -resume
2016-07-25 09:58:01\t5 min\tmodest_bartik\tERR\t6b9515aba6\t5910a50f-8656-4765-aa79-f07cef912062\t.nextflow run hello
'''

        when:
        when = HistoryFile.TIMESTAMP_FMT.parse('2016-07-24 17:43:34')
        history.update('small_cirum',false,when)
        then:
        history.text == '''
2016-07-24 16:43:16\t10m\tevil_pike\tOK\t6b9515aba6\te710da1b-ce06-482f-bbcf-987a507f85d1\t.nextflow run hello
2016-07-24 16:43:34\t-\tgigantic_keller\t-\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello
2016-07-24 16:43:34\t1h\tsmall_cirum\tERR\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello -resume
2016-07-25 09:58:01\t5 min\tmodest_bartik\tERR\t6b9515aba6\t5910a50f-8656-4765-aa79-f07cef912062\t.nextflow run hello
'''

        when:
        when = HistoryFile.TIMESTAMP_FMT.parse('2016-07-24 16:43:50')
        history.update('gigantic_keller',true,when)
        then:
        history.text == '''
2016-07-24 16:43:16\t10m\tevil_pike\tOK\t6b9515aba6\te710da1b-ce06-482f-bbcf-987a507f85d1\t.nextflow run hello
2016-07-24 16:43:34\t16s\tgigantic_keller\tOK\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello
2016-07-24 16:43:34\t1h\tsmall_cirum\tERR\t6b9515aba6\t5a6d3877-8823-4ed6-b7fe-2b6748ed4ff9\t.nextflow run hello -resume
2016-07-25 09:58:01\t5 min\tmodest_bartik\tERR\t6b9515aba6\t5910a50f-8656-4765-aa79-f07cef912062\t.nextflow run hello
'''

    }

}
