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

import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MailTest extends Specification {

    def 'should capture mail params' () {

        given:
        def closure = {
            from 'jim@dot.com'
            to 'paolo@dot.com'
            cc 'you@dot.com'
            bcc 'mrhide@dot.com'
            type 'text/html'
            subject 'Hi there'
            charset 'utf-8'
            body 'Hello world'
            attach 'foo.png'
            attach (['this.txt','that.txt'])

        }

        when:
        def mail = new Mail()
        mail.with(closure)
        then:
        mail.from == 'jim@dot.com'
        mail.to == 'paolo@dot.com'
        mail.cc == 'you@dot.com'
        mail.bcc == 'mrhide@dot.com'
        mail.type == 'text/html'
        mail.subject ==  'Hi there'
        mail.charset == 'utf-8'
        mail.body == 'Hello world'
        mail.attachments == [new Attachment('foo.png'), new Attachment('this.txt'), new Attachment('that.txt')]

    }


    def 'should add attachments' () {
        given:
        Mail mail

        when:
        mail = new Mail()
        mail.attach('/some/file.txt')
        then:
        mail.attachments == [new Attachment('/some/file.txt')]

        when:
        mail = new Mail()
        mail.attach(new File('x.txt'))
        then:
        mail.attachments == [new Attachment('x.txt')]

        when:
        mail = new Mail()
        mail.attach(Paths.get('x.txt'))
        then:
        mail.attachments == [new Attachment('x.txt')]

        when:
        mail = new Mail()
        mail.attach("file.${1}")
        then:
        mail.attachments == [new Attachment('file.1')]

        when:
        mail = new Mail()
        mail.attach(['foo.txt','bar.txt'])
        then:
        mail.attachments == [new Attachment('foo.txt'), new Attachment('bar.txt')]

        when:
        mail = new Mail()
        mail.attach 'pic.png', contentId: 'my-img'
        then:
        mail.attachments == [new Attachment('pic.png', contentId:'my-img')]

        when:
        mail = new Mail()
        mail.attach( new Attachment('/file.txt') )
        then:
        mail.attachments ==  [new Attachment('/file.txt')]
        when:
        mail = new Mail()
        mail.attach(new Object())
        then:
        thrown(IllegalArgumentException)
    }

    def 'should create a mail from a Map' () {

        given:
        def map = [
                from:'me@google.com',
                to: 'you@nextflow.com',
                cc: 'hola@dot.com, hello@dot.com',
                bcc: 'foo@host.com',
                subject: 'this is a notification',
                charset: 'utf-8',
                type: 'text/html',
                body: 'Hello world',
                text: 'Pura vida',
                attach: '/some/file'
        ]

        when:
        def mail = Mail.of(map)
        then:
        mail.from == 'me@google.com'
        mail.to == 'you@nextflow.com'
        mail.cc == 'hola@dot.com, hello@dot.com'
        mail.bcc == 'foo@host.com'
        mail.subject == 'this is a notification'
        mail.charset == 'utf-8'
        mail.type == 'text/html'
        mail.body == 'Hello world'
        mail.text == 'Pura vida'
        mail.attachments == [new Attachment('/some/file')]
    }
}
