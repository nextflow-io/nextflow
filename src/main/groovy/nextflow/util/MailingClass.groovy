package nextflow.util


import groovy.util.logging.Slf4j

import javax.mail.*;
import javax.mail.internet.*;


@Slf4j
class MailingClass {

    // Get system properties
    private final Properties properties = System.getProperties();
    private final String imapHost;
    private final int imapPort;

    private final String smtpHost;
    String osName = System.getProperty("os.name").toLowerCase();

    protected MailingClass(String smtpHost, int smtpPort,
                        String imapHost, int imapPort,
                        String user, String password) {

        this.smtpHost = smtpHost;
        this.imapHost = imapHost;
        this.imapPort = imapPort;
        properties.setProperty("mail.smtp.host", smtpHost);     //host
        properties.setProperty("mail.smtp.port", Integer.toString(smtpPort));
        properties.setProperty("mail.user", user);
        properties.setProperty("mail.password", password);
        properties.setProperty("mail.store.protocol", "imap");

    }
    protected Message createJavaMessage(String to, String from, String subject, String body, String attachment) {

        // Get the default Session object.
        def session = javax.mail.Session.getDefaultInstance(properties);

        // Create a default MimeMessage object.
        MimeMessage message = new MimeMessage(session);

        // Set From: header field of the header.
        message.setFrom(new InternetAddress(from));

        // Set To: header field of the header.
        message.addRecipient(Message.RecipientType.TO, new InternetAddress(to));

        // Set Subject: header field
        message.setSubject(subject);

        // Create the message part
        BodyPart messageBodyPart = new MimeBodyPart();

        // Fill the message
        messageBodyPart.setText(body);

        // Create a multipar message
        Multipart multipart = new MimeMultipart();

        // Set text message part
        multipart.addBodyPart(messageBodyPart);

        if(attachment!=null)
        {
            // Part two is attachment
            messageBodyPart = new MimeBodyPart();
            messageBodyPart.attachFile(new File(attachment));
            multipart.addBodyPart(messageBodyPart);
        }
        // Send the complete message parts
        message.setContent(multipart);
        return message
    }

    protected void sendByJavaMail(String to, String from, String subject, String body, String attachment) {

        def message = createJavaMessage(to, from, subject, body, attachment)
        // Send message
        Transport.send(message);
    }

    protected List<String> sendBySysMail(String to, String from, String subject, String body, String attachment) {

        List<String> result = new LinkedList<String>()

        if (osName=="mac" ) {
            //command="[\'sendmail\',\'-a\']"
            result.add("sendmail")
            result.add(to)
            result.add(body)
        }else {
            //command=["mail", "-s", subject, to]
            result.add("mail")
            result.add("-s")
            result.add(subject)
            result.add(to)
            result.add(" <<< "+body)
            // MAIL FROM A FILE
            //  mail -s "This is Subject" someone@example.com < /path/to/file
            // ATTACHMENT
            //  "This is message body" | mail -s "This is Subject" -r "Harry<harry@gmail.com>" -a /path/to/file someone@example.com
        }
        return result;
    }

    public void send(String to, String from, String subject, String body, String attachment) {
        if( this.smtpHost == 'localhost' ) {
            sendByJavaMail(to, from, subject, body, attachment)
        }
        else {
            def command = sendBySysMail(to, from, subject, body, attachment)

            def process = new ProcessBuilder().command(command).start()

            def result = process.waitFor()

            if( result != 0 ) {
                log.error "Unable to send mail -- command executed: $command"
            }
        }
    }

    public Message[] receiveMail(String user, String password) throws MessagingException {
        javax.mail.Session session = javax.mail.Session.getDefaultInstance(properties);

        Store store = session.getStore("imap");
        store.connect(imapHost, imapPort, user, password);
        Folder inbox = store.getFolder("INBOX");

        inbox.open(Folder.READ_ONLY);

        return inbox.getMessages();
    }
}