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
        mail.attachments == [new File('foo.png'), new File('this.txt'), new File('that.txt')]

    }


    def 'should add attachments' () {
        given:
        Mail mail

        when:
        mail = new Mail()
        mail.attach('/some/file.txt')
        then:
        mail.attachments == [new File('/some/file.txt')]

        when:
        mail = new Mail()
        mail.attach(new File('x.txt'))
        then:
        mail.attachments == [new File('x.txt')]

        when:
        mail = new Mail()
        mail.attach(Paths.get('x.txt'))
        then:
        mail.attachments == [new File('x.txt')]

        when:
        mail = new Mail()
        mail.attach("file.${1}")
        then:
        mail.attachments == [new File('file.1')]

        when:
        mail = new Mail()
        mail.attach(['foo.txt','bar.txt'])
        then:
        mail.attachments == [new File('foo.txt'), new File('bar.txt')]

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
        mail.attachments == [new File('/some/file')]
    }
}
