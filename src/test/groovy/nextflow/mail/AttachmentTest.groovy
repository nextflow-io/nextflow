/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
