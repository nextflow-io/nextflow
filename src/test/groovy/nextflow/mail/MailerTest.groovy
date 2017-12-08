package nextflow.mail
import javax.mail.Message
import javax.mail.internet.InternetAddress
import javax.mail.internet.MimeMessage
import javax.mail.internet.MimeMultipart
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import com.icegreen.greenmail.util.GreenMail
import com.icegreen.greenmail.util.ServerSetupTest
import spock.lang.Specification

class MailerTest extends Specification {


    def 'should return config properties'() {
        when:
        def SMTP = [host: 'google.com', port: '808', user: 'foo', password: 'bar']
        def mailer = new Mailer( config: [smtp: SMTP, other: 1]  )
        def props = mailer.createProps()

        then:
        props.get('mail.smtp.user') == 'foo'
        props.get('mail.smtp.password') == 'bar'
        props.get('mail.smtp.host') == 'google.com'
        props.get('mail.smtp.port') == '808'
        !props.containsKey('mail.other')

    }


    def "sending mails using javamail"() {

        given:
        final USER = 'foo'
        final PASSWORD = 'secret'
        final EMAIL = 'yo@nextflow.com'
        final green = new GreenMail(ServerSetupTest.SMTP)
        green.setUser(EMAIL, USER, PASSWORD)
        green.start()

        def SMTP = [ host: '127.0.0.1', port: green.smtp.port, user: USER, password: PASSWORD]
        def mailer = new Mailer( config: [smtp: SMTP])

        String TO = "receiver@testmailingclass.net"
        String FROM = 'paolo@gmail.com'
        String SUBJECT = "Sending test"
        String CONTENT = "This content should be sent by the user."

        when:
        mailer.to = TO
        mailer.from = FROM
        mailer.subject = SUBJECT
        mailer.content = CONTENT
        mailer.send()

        then:
        green.receivedMessages.size() == 1
        Message message = green.receivedMessages[0]
        message.from == [new InternetAddress(FROM)]
        message.allRecipients.contains(new InternetAddress(TO))
        message.subject == SUBJECT
        message.getContent() instanceof MimeMultipart
        (message.getContent() as MimeMultipart).getBodyPart(0).content == CONTENT

        cleanup:
        green?.stop()
    }

    def "sending mails using java with attachment"() {

        given:
        final USER = 'foo'
        final PASSWORD = 'secret'
        final EMAIL = 'yo@nextflow.com'
        final green = new GreenMail(ServerSetupTest.SMTP)
        green.setUser(EMAIL, USER, PASSWORD)
        green.start()

        def SMTP = [ host: '127.0.0.1', port: green.smtp.port, user: USER, password: PASSWORD]
        def mailer = new Mailer( config: [smtp: SMTP])

        String TO = "receiver@testmailingclass.net"
        String FROM = 'paolo@nextflow.io'
        String SUBJECT = "Sending test"
        String CONTENT = "This content should be sent by the user."
        Path ATTACH = Files.createTempFile('test', null)
        ATTACH.text = 'This is the file attachment content'

        when:
        mailer.from = FROM
        mailer.to = TO
        mailer.subject = SUBJECT
        mailer.content = CONTENT
        mailer.addAttachment( ATTACH )
        mailer.send()

        then:
        green.receivedMessages.size() == 1
        Message message = green.receivedMessages[0]
        message.from == [new InternetAddress(FROM)]
        message.allRecipients.contains(new InternetAddress(TO))
        message.subject == SUBJECT
        (message.getContent() as MimeMultipart).getCount() == 2
        //(message.getContent() as MimeMultipart).getBodyPart(1).getContent() == ''

        cleanup:
        ATTACH?.delete()
        green?.stop()
    }


    def 'should send with java' () {

        given:
        def mailer = Spy(Mailer)
        def MSG = Mock(MimeMessage)

        when:
        mailer.config = [smtp: [host:'foo.com'] ]
        mailer.send()
        then:
        1 * mailer.createMimeMessage() >> MSG
        1 * mailer.sendViaJavaMail(MSG) >> null

    }

    def 'should send with command line' () {
        given:
        def mailer = Spy(Mailer)
        def MSG = Mock(MimeMessage)

        when:
        mailer.send()
        then:
        1 * mailer.createMimeMessage() >> MSG
        1 * mailer.sendViaSendmail(MSG) >> null
    }


    def 'should create mime message' () {

        given:
        MimeMessage msg

        when:
        msg = new Mailer(from:'foo@gmail.com').createMimeMessage()
        then:
        msg.getFrom().size()==1
        msg.getFrom()[0].toString() == 'foo@gmail.com'

        when:
        msg = new Mailer(from:'one@gmail.com, two@google.com').createMimeMessage()
        then:
        msg.getFrom().size()==2
        msg.getFrom()[0].toString() == 'one@gmail.com'
        msg.getFrom()[1].toString() == 'two@google.com'

        when:
        msg = new Mailer(to:'foo@gmail.com, bar@google.com').createMimeMessage()
        then:
        msg.getRecipients(Message.RecipientType.TO).size()==2
        msg.getRecipients(Message.RecipientType.TO)[0].toString() == 'foo@gmail.com'
        msg.getRecipients(Message.RecipientType.TO)[1].toString() == 'bar@google.com'

        when:
        msg = new Mailer(cc:'foo@gmail.com, bar@google.com').createMimeMessage()
        then:
        msg.getRecipients(Message.RecipientType.CC).size()==2
        msg.getRecipients(Message.RecipientType.CC)[0].toString() == 'foo@gmail.com'
        msg.getRecipients(Message.RecipientType.CC)[1].toString() == 'bar@google.com'

        when:
        msg = new Mailer(bcc:'one@gmail.com, two@google.com').createMimeMessage()
        then:
        msg.getRecipients(Message.RecipientType.BCC).size()==2
        msg.getRecipients(Message.RecipientType.BCC)[0].toString() == 'one@gmail.com'
        msg.getRecipients(Message.RecipientType.BCC)[1].toString() == 'two@google.com'

        when:
        msg = new Mailer(subject: 'this is a test', content: 'Ciao mondo').createMimeMessage()
        then:
        msg.getSubject() == 'this is a test'
        msg.getContent() instanceof MimeMultipart
        msg.getContent().getCount() == 1
        msg.getContent().getBodyPart(0).getContent() == 'Ciao mondo'
    }

    def 'should add attachments' () {
        given:
        Mailer mail

        when:
        mail = new Mailer().addAttachment('/some/file.txt')
        then:
        mail.attachments == [new File('/some/file.txt')]

        when:
        mail = new Mailer().addAttachment(new File('x.txt'))
        then:
        mail.attachments == [new File('x.txt')]

        when:
        mail = new Mailer().addAttachment(Paths.get('x.txt'))
        then:
        mail.attachments == [new File('x.txt')]

        when:
        mail = new Mailer().addAttachment("file.${1}")
        then:
        mail.attachments == [new File('file.1')]

        when:
        mail = new Mailer().addAttachment(['foo.txt','bar.txt'])
        then:
        mail.attachments == [new File('foo.txt'), new File('bar.txt')]

        when:
        new Mailer().addAttachment(new Object())
        then:
        thrown(IllegalArgumentException)
    }

    def 'should fetch config properties' () {

        given:
        def ENV = [NXF_SMTP_USER: 'jim', NXF_SMTP_PASSWORD: 'secret', NXF_SMTP_HOST: 'g.com', NXF_SMTP_PORT: '864']
        def SMTP = [host:'hola.com', user:'foo', password: 'bar', port: 234]
        Mailer mail

        when:
        mail = new Mailer(config: [smtp: SMTP])
        then:
        mail.host == 'hola.com'
        mail.user == 'foo'
        mail.password == 'bar'
        mail.port == 234

        when:
        mail = new Mailer(config: [smtp: [host: 'local', port: '999']], env: ENV)
        then:
        mail.host == 'local'
        mail.port == 999
        mail.user == 'jim'
        mail.password == 'secret'

        when:
        mail = new Mailer(env: ENV)
        then:
        mail.host == 'g.com'
        mail.port == 864
        mail.user == 'jim'
        mail.password == 'secret'
    }

    def 'should config the mailer' () {

        given:
        def CFG = '''
            mail {
              from = 'paolo@nf.com'
              subject = 'nextflow notification'

              smtp.host = 'foo.com'
              smtp.port = '43'
              smtp.user = 'jim'
              smtp.password = 'yo'
            }
            '''

        when:
        def mailer = new Mailer()
        mailer.config = new ConfigSlurper().parse(CFG).mail
        then:
        mailer.getUser() == 'jim'
        mailer.getPassword() == 'yo'
        mailer.getPort() == 43
        mailer.getHost() == 'foo.com'

        mailer.from == 'paolo@nf.com'
        mailer.subject == 'nextflow notification'
    }


    def 'should capture send params' () {

        given:
        def mailer = Spy(Mailer)

        when:
        mailer.send {
            to 'paolo@dot.com'
            from 'yo@dot.com'
            subject 'This is a test'
            content 'Hello there'
        }

        then:
        1 * mailer.send( [to: 'paolo@dot.com', from:'yo@dot.com', subject: 'This is a test', content: 'Hello there'] ) >> null

    }


}
