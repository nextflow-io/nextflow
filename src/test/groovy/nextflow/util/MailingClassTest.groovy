package nextflow.util

import com.dumbster.smtp.SimpleSmtpServer
import com.dumbster.smtp.SmtpMessage
import spock.lang.Specification

import com.icegreen.greenmail.user.GreenMailUser
import com.icegreen.greenmail.util.GreenMail
import com.icegreen.greenmail.util.ServerSetupTest

import javax.mail.Message
import javax.mail.internet.InternetAddress
import javax.mail.internet.MimeMessage

/**
 * Created by edgar on 8/11/17.
 */
class MailingClassTest extends Specification{

    GreenMail greenMail = new GreenMail(ServerSetupTest.SMTP_IMAP)
    GreenMailUser greenMailUser = greenMail.setUser("someuser@somewhere.com", "someUser", "somePassword")

    def setup() {
        greenMail.start()
    }

    def cleanup() {
        greenMail.stop()
    }

    def "sending mails using JAVA"() {

        given:
        String smtpHost = 'localhost'

        MailingClass mailingClass = new MailingClass(
                smtpHost, ServerSetupTest.SMTP.port,
                "localhost", ServerSetupTest.IMAP.port,
                greenMailUser.login, greenMailUser.password)

        String to = "receiver@testmailingclass.net"
        String from = greenMailUser.email
        String subject = "Sending test"
        String content = "This content should be sent by the user."
        String attachment = "/home/edgar/Desktop/nextflowGit/src/test/groovy/nextflow/util/example.jpg"

        when:
        mailingClass.send(to, from, subject, content, attachment)

        then:
        greenMail.receivedMessages.size() == 1
        Message message = greenMail.receivedMessages[0]
        message.from == [new InternetAddress(from)]
        message.allRecipients.contains(new InternetAddress(to))
        message.subject == subject
        //message.content == "${content}\r\n"
    }
    def "sending mails using LINUX/MAC"() {

        given:
        String smtpHost = ''

        MailingClass mailingClass = new MailingClass(
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

        MailingClass mailingClass = new MailingClass(
                smtpHost, ServerSetupTest.SMTP.port,
                "localhost", ServerSetupTest.IMAP.port,
                greenMailUser.login, greenMailUser.password)

        String from = "sender@testmailingclass.net"
        String subject = "Sending test"
        String content = "This content should be received by the user."
        deliverMessage(from, subject, content)

        when:
        Message[] messages = mailingClass.receiveMail(
                greenMailUser.login, greenMailUser.password)

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
