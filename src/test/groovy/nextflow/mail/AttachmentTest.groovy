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

package nextflow.mail

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AttachmentTest extends Specification {


    def 'should create resource attachment' () {

        when:
        def attach = Attachment.resource('foo/bar', contentId: 'the-cid')
        then:
        attach.file == null
        attach.resource == 'foo/bar'
        attach.contentId == 'the-cid'

    }

    def 'should crate attachment'  () {

        given:
        Attachment attach

        when:
        attach = new Attachment('/some/path/foo.png')
        then:
        attach.file == new File('/some/path/foo.png')
        attach.fileName == 'foo.png'
        attach.contentId == null
        attach.description == null
        attach.disposition == null

        when:
        attach = new Attachment('/some/path/foo.png', contentId: 'id-1', description: 'Hola', fileName: 'bar.png', disposition: 'inline')
        then:
        attach.file == new File('/some/path/foo.png')
        attach.fileName == 'bar.png'
        attach.description == 'Hola'
        attach.contentId == 'id-1'
        attach.disposition == 'inline'

        when:
        attach = Attachment.resource('jar:/some/path/foo.png', contentId: '<foo>')
        then:
        attach.file == null
        attach.resource == 'jar:/some/path/foo.png'
        attach.fileName == 'foo.png'
        attach.contentId == '<foo>'

        when:
        attach = Attachment.resource([:], 'jar:foo.png')
        then:
        attach.file == null
        attach.resource == 'jar:foo.png'
        attach.fileName == 'foo.png'

    }

}
