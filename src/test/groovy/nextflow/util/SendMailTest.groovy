package nextflow.util
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
/**
 * Created by edgar on 8/11/17.
 */
class SendMailTest extends Specification{

    final greenMail = new GreenMail(ServerSetupTest.SMTP_IMAP)
    final greenMailUser = greenMail.setUser("someuser@somewhere.com", "someUser", "somePassword")

    def setup() {
        greenMail.start()
    }

    def cleanup() {
        greenMail.stop()
    }

    def 'should return config properties' () {

        when:
        def mailer = new SendMail(host: 'google.com', protocol: 'imap', port: '808', user: 'foo', password: 'bar')
        def props = mailer.createProps()
        then:
        props.get('mail.user') == 'foo'
        props.get('mail.password') == 'bar'
        props.get('mail.imap.host') == 'google.com'
        props.get('mail.imap.port') == '808'
        props.get('mail.store.protocol') == 'imap'

        when:
        mailer = new SendMail(host: 'crg.com', protocol: 'smtp', port: '25', user: 'xxx', password: 'zzz')
        props = mailer.createProps()
        then:
        props.get('mail.user') == 'xxx'
        props.get('mail.password') == 'zzz'
        props.get('mail.smtp.host') == 'crg.com'
        props.get('mail.smtp.port') == '25'
        props.get('mail.store.protocol') == 'smtp'

        when:
        mailer = new SendMail(protocol: 'foo')
        mailer.createProps()
        then:
        thrown(IllegalArgumentException)
    }

    def 'should get send mail command line' () {

        given:
        def ATTACHMENT = Paths.get('/some/attachment.txt')

        when:
        def mailer = new SendMail(mailer:'mail')
        def cli = mailer.sendBySysMail('paolo@crg.eu', 'edgar@crg.eu', 'hello', 'mail content', ATTACHMENT)
        then:
        cli == []

        when:
        mailer = new SendMail(mailer:'sendmail')
        cli = mailer.sendBySysMail('paolo@crg.eu', 'edgar@crg.eu', 'hello', 'mail content', ATTACHMENT)
        then:
        cli == []

    }

    def "sending mails using JAVA"() {

        given:
        def mailingClass = new SendMail(
                host: 'localhost',
                port: ServerSetupTest.SMTP.port as String,
                user: greenMailUser.login,
                password: greenMailUser.password)

        String TO = "receiver@testmailingclass.net"
        String FROM = greenMailUser.email
        String SUBJECT = "Sending test"
        String CONTENT = "This content should be sent by the user."

        when:
        mailingClass.send(TO, FROM, SUBJECT, CONTENT)

        then:
        greenMail.receivedMessages.size() == 1
        Message message = greenMail.receivedMessages[0]
        message.from == [new InternetAddress(FROM)]
        message.allRecipients.contains(new InternetAddress(TO))
        message.subject == SUBJECT
        println message.content.class
        (message.content as MimeMultipart).getBodyPart(1) == ''
    }

    def "sending mails using JAVA with attachment"() {

        given:
        def temp = Files.createTempFile('test',null)
        temp.text = 'attachment content blah blah'

        def mailingClass = new SendMail(
                host: 'localhost',
                port: ServerSetupTest.SMTP.port as String,
                user: greenMailUser.login,
                password: greenMailUser.password)

        String TO = "receiver@testmailingclass.net"
        String FROM = greenMailUser.email
        String SUBJECT = "Sending test"
        String CONTENT = "This content should be sent by the user."
        Path ATTACHMENT = temp

        when:
        mailingClass.send(TO, FROM, SUBJECT, CONTENT)

        then:
        greenMail.receivedMessages.size() == 1
        Message message = greenMail.receivedMessages[0]
        message.from == [new InternetAddress(FROM)]
        message.allRecipients.contains(new InternetAddress(TO))
        message.subject == SUBJECT
        println message.content.class
        (message.content as MimeMultipart).getBodyPart(1) == ''


        cleanup:
        temp?.delete()
    }


    def "sending mails using LINUX/MAC"() {

        given:
        String smtpHost = ''

        SendMail mailingClass = new SendMail(
                smtpHost, ServerSetupTest.SMTP.port,
                "localhost", ServerSetupTest.IMAP.port,
                greenMailUser.login, greenMailUser.password)

        String to = "receiver@testmailingclass.net"
        String from = greenMailUser.email
        String subject = "Sending test"
        String content = "This content should be sent by the user."
        String attachment = "/home/edgar/Desktop/nextflowGit/src/test/groovy/nextflow/util/example.jpg"

        when:
        mailingClass.osName='mac'
        mailingClass.send(to, from, subject, content, attachment)

        then:
        greenMail.receivedMessages.size() == 1
        Message message = greenMail.receivedMessages[0]
        message.from == [new InternetAddress(from)]
        message.allRecipients.contains(new InternetAddress(to))
        message.subject == subject

    }
    def "receiving mails using JAVA"() {
        given:
        String smtpHost = 'localhost'

        SendMail mailingClass = new SendMail(
                smtpHost, ServerSetupTest.SMTP.port,
                "localhost", ServerSetupTest.IMAP.port,
                greenMailUser.login, greenMailUser.password)

        String from = "sender@testmailingclass.net"
        String subject = "Sending test"
        String content = "This content should be received by the user."
        deliverMessage(from, subject, content)

        when:
        Message[] messages = mailingClass.receiveMail(greenMailUser.login, greenMailUser.password)

        then:
        messages.size() == 1
        Message message = messages[0]
        message.from == [new InternetAddress(from)]
        message.allRecipients == [new InternetAddress(greenMailUser.email)]
        message.subject == subject


    }


    private MimeMessage deliverMessage(String from, String subject, String text) {
        MimeMessage message = new MimeMessage((javax.mail.Session) null)
        message.setFrom(new InternetAddress(from))
        message.setRecipients(Message.RecipientType.TO, greenMailUser.email)
        message.setSubject(subject)
        message.setText(text)
        greenMailUser.deliver(message)
        return message
    }

}
