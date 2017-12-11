package nextflow.util

import javax.mail.Header
import javax.mail.Message
import javax.mail.Part
import javax.mail.internet.InternetAddress
import javax.mail.internet.MimeBodyPart
import javax.mail.internet.MimeMultipart
import java.nio.file.Path
import com.icegreen.greenmail.util.GreenMail
import com.icegreen.greenmail.util.ServerSetupTest
import spock.lang.Specification

import java.nio.file.Paths

/**
 * Created by edgar on 8/11/17.
 */
class SendMailTest extends Specification {

    final greenMail = new GreenMail(ServerSetupTest.SMTP_IMAP)
    final greenMailUser = greenMail.setUser("someuser@somewhere.com", "someUser", "somePassword")

    def setup() {
        greenMail.start()
    }

    def cleanup() {
        greenMail.stop()
    }

    def 'should return config properties'() {
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

    def 'should get send mail command line'() {
        given:
        def body = 'mail content'
        def subject = 'subject'
        def to = 'paolo@crg.eu'
        def from = 'edgar@crg.es'
        def ATTACHMENT = Paths.get('/home/edgar/Desktop/nextflow/src/test/groovy/nextflow/util/email.html')
        List<String> resultMail = ['echo', body, '|', 'mail', '-s', subject, '-a', ATTACHMENT, to]
        def auxBody = body+"\n ATTACHMENT FILE: \n"+ATTACHMENT.text
        List<String> resultMailNoFile = ['echo', auxBody, '|', 'mail', '-s', subject, to]
        List<String> resultSendMail = ['sendmail', '-t', to, '<',ATTACHMENT]

        when:
        def mailer = new SendMail(mailer: 'mail', includeFile:true)
        def cli = mailer.sendBySysMail(to, from, subject, body, ATTACHMENT)
        then:
        cli == resultMail

        when:
        mailer = new SendMail(mailer: 'mail', includeFile:false)
        cli = mailer.sendBySysMail(to, from, subject, body, ATTACHMENT)
        then:
        cli == resultMailNoFile

        when:
        mailer = new SendMail(mailer: 'sendmail', includeFile:true)
        cli = mailer.sendBySysMail(to, from, subject, body, ATTACHMENT)
        then:
        cli == resultSendMail

        when:
        mailer = new SendMail(mailer: 'sendmail', includeFile:false)
        cli = mailer.sendBySysMail(to, from, subject, body, ATTACHMENT)
        then:
        cli == resultSendMail
    }

    def "sending mails using JAVA without ATTACHMENT"() {

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
    }

    def "sending mails using JAVA with ATTACHMENT W/O file"() {
        given:
        def attachmentList = new ArrayList<String>()

        def mailingClass = new SendMail(
                host: 'localhost',
                port: ServerSetupTest.SMTP.port as String,
                user: greenMailUser.login,
                password: greenMailUser.password)
        mailingClass.includeFile=true
        String TO = "receiver@testmailingclass.net"
        String FROM = greenMailUser.email
        String SUBJECT = "Sending test"
        String CONTENT = "This content should be sent by the user."
        Path ATTACHMENT = Paths.get("/home/edgar/Desktop/hello.html")
        String fileAttached = new File(ATTACHMENT.toString()).text
        attachmentList.add(fileAttached)

        when:
        mailingClass.send(TO, FROM, SUBJECT, CONTENT, ATTACHMENT)

        then:
        greenMail.receivedMessages.size() == 1
        Message message = greenMail.receivedMessages[0]
        message.from == [new InternetAddress(FROM)]
        message.allRecipients.contains(new InternetAddress(TO))
        message.subject == SUBJECT
        //Need replaceAll to remove \r and make equal the comparison.
        this.getAttachmentAsString(message).get(0).toString().replaceAll("\\r","") == attachmentList.get(0).toString()

        greenMail.reset()

        when:
        mailingClass.includeFile=false
        mailingClass.send(TO, FROM, SUBJECT, CONTENT, ATTACHMENT)

        then:
        greenMail.receivedMessages.size() == 1
        Message message2 = greenMail.receivedMessages[0]
        message2.from == [new InternetAddress(FROM)]
        message2.allRecipients.contains(new InternetAddress(TO))
        message2.subject == SUBJECT

        this.getBodyAsString(message2).get(1).replaceAll("\\r","") == attachmentList.get(0).toString()
    }

    private void printMessage(Message message) {
        if (message.getContent() instanceof MimeMultipart) {
            System.out.println("-------------------------------------------------------");
            MimeMultipart mimeMultipart = (MimeMultipart) message.getContent();
            Enumeration headers = message.getAllHeaders();
            println("*** Headers ***")
            while (headers.hasMoreElements()) {
                Header h = (Header) headers.nextElement();
                System.out.println(h.getName() + " : " + h.getValue());
            }
            println("*** Headers END***")
            for (int i = 0; i < mimeMultipart.getCount(); i++) {
                MimeBodyPart part = (MimeBodyPart) mimeMultipart.getBodyPart(i);
                if (Part.ATTACHMENT.equalsIgnoreCase(part.getDisposition())) {
                    println("***Attachment ***")
                    println("File name: " + part.getFileName())

                    println(part.getContent().toString())
                    println("***Attachment END***")
                } else {
                    println("Content: ")
                    System.out.println(mimeMultipart.getBodyPart(i).getContentType() + " : " + mimeMultipart.getBodyPart(i).getContent());
                }
            }
        } else {
            System.out.println("else " + message.getContent());
        }
        System.out.println("-------------------------------------------------------");
    }

    private ArrayList<String> getBodyAsString(Message message) {
        List<String> result = new ArrayList<String>();

        MimeMultipart mimeMultipart = (MimeMultipart) message.getContent();

        for (int i = 0; i < mimeMultipart.getCount(); i++) {
            MimeBodyPart part = (MimeBodyPart) mimeMultipart.getBodyPart(i);
            if(!Part.ATTACHMENT.equalsIgnoreCase(part.getDisposition())){
                result.add(part.getContent().toString())
            }
        }
        return result;
    }

    private ArrayList<String> getAttachmentAsString(Message message) {
        List<String> result = new ArrayList<String>();

        MimeMultipart mimeMultipart = (MimeMultipart) message.getContent();

        for (int i = 0; i < mimeMultipart.getCount(); i++) {
            MimeBodyPart part = (MimeBodyPart) mimeMultipart.getBodyPart(i);
            if (Part.ATTACHMENT.equalsIgnoreCase(part.getDisposition())) {
                part.toString();
                result.add(part.getContent().toString())
            }
        }
        return result;
    }

}