/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.util

import java.nio.file.Path
import javax.mail.*;
import javax.mail.internet.*;

import groovy.transform.Memoized
import groovy.util.logging.Slf4j

/**
 * This class implement the send mail functionality
 *
 * @author Edgar Garriga
 */
@Slf4j
class SendMail {

    static final List<String> PROTOCOLS = ['smtp','imap']

    /** the type of the server ie. IMAP or SMTP */
    private String protocol = 'smtp'

    private String host;

    private String port;

    private String user

    private String password

    /**To send the attachment as a file attached or as a text in the body **/
    private boolean includeFile

    /** Either system tool `sendmail` or `mail` */
    private String mailer

    private String SYS_MAILER

    {
        SYS_MAILER = findSysMailer()
    }

    @Memoized
    protected String getMailer() {
        mailer ?: findSysMailer()
    }

    protected String findSysMailer() {

        // first try `sendmail`
        if( runCommand("command -v sendmail &>/dev/null") == 0 )
            return 'sendmail'

        else if( runCommand("command -v mail &>/dev/null") == 0  )
            return 'mail'

        log.warn "Missing system mail command -- Make sure to have either `sendmail` or `mail` tool installed otherwise configure your mail server properties in the nextflow config file"
        return null
    }

    protected int runCommand(String cmd) {
        def proc = ['bash','-c',cmd].execute()
        proc.waitForOrKill(1_000)
        return proc.exitValue()
    }

    protected Properties createProps() {
        def properties = System.getProperties();

        def type = protocol.toLowerCase()
        if( !PROTOCOLS.contains(type) )
            throw new IllegalArgumentException("Not a valid email protocol: $type")

        properties.setProperty("mail.${type}.host", host);     //host
        properties.setProperty("mail.${type}.port", port);
        properties.setProperty("mail.store.protocol", type);

        properties.setProperty("mail.user", user);
        properties.setProperty("mail.password", password);

        return properties
    }

    protected Message createJavaMessage(String to, String from, String subject, String body, Path attachment=null) {

        // Get the default Session object.
        def session = javax.mail.Session.getDefaultInstance(createProps());

        // Create a default MimeMessage object.
        def message = new MimeMessage(session);

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

        // Create a multipart message
        Multipart multipart = new MimeMultipart();

        // Set text message part
        multipart.addBodyPart(messageBodyPart);

        if(attachment && includeFile)
        {
            //println("*****FILE IN ATTACHMENT")
            // Part two is attachment
            messageBodyPart = new MimeBodyPart();
            messageBodyPart.attachFile(new File(attachment.toString()));
            multipart.addBodyPart(messageBodyPart);
        } else if(attachment && !includeFile){          //send attachment in body
            //println("*****NOT ***FILE IN ATTACHMENT")

            messageBodyPart = new MimeBodyPart();
            String fileString = new File(attachment.toString()).text
            //println("fileString: \n"+fileString+"+++++")
            messageBodyPart.setText(fileString)
            multipart.addBodyPart(messageBodyPart)
        }
        // Send the complete message parts
        message.setContent(multipart);
        return message
    }

    protected void sendByJavaMail(String to, String from, String subject, String body, Path attachment=null) {

        def message = createJavaMessage(to, from, subject, body, attachment)
        // Send message
        Transport.send(message);
    }

    protected List<String> sendBySysMail(String to, String from, String subject, String body, Path attachment=null) {

        List<String> result = new LinkedList<String>()

        def mailer = getMailer()

        if ( mailer=="sendmail" ) {
            //  /usr/sbin/sendmail -t edgar.garriga@crg.es < /users/cn/egarriga/example.txt
            /*
            Subject: Subject INSIDE
            From: a@example.com
            ---
            THIS IS AN ATTACHMENT
            hola mundo
            attachment
            ---
            */
            result.add("sendmail")
            result.add("-t")
            result.add(to)
            result.add("<")
            result.add(attachment)
        }
        else if( mailer=='mail' ) {
            if(includeFile){
                //echo "BODY" | mail -s "SUBJECT"  -a /users/cn/egarriga/example.html edgar.garriga@crg.es
                result.add("echo")
                result.add(body)
                result.add("|")
                result.add("mail")
                result.add("-s")
                result.add(subject)
                result.add("-a")
                result.add(attachment)
                result.add(to)
            }else{
                //echo "BODY" | mail -s "SUBJECT" edgar.garriga@crg.es
                String attachmentString = new File(attachment.toString()).text
                body = body+"\n ATTACHMENT FILE: \n"+attachmentString

                result.add("echo")
                result.add(body)
                result.add("|")
                result.add("mail")
                result.add("-s")
                result.add(subject)
                result.add(to)
            }
        }
            //todo check MUTT

            //todo sendEmail
        //http://caspian.dotconf.net/menu/Software/SendEmail/
        /*  result.add("sendEmail")
            result.add("-f")
            result.add(from)
            result.add("-t")
            result.add(to)
            result.add("-u")
            result.add(subject)
            result.add("-m")
            result.add(body)
            result.add("-a")
            result.add(attachment)
        */
        else
            throw new IllegalStateException("Invalid mailer function: $mailer")

        return result;
    }

    public void send(String to, String from, String subject, String body, Path attachment=null) {
        if(!attachment)
            includeFile=false
        if( host ) {
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
        javax.mail.Session session = javax.mail.Session.getDefaultInstance(createProps());

        Store store = session.getStore("imap");
        store.connect(host, port, user, password);
        Folder inbox = store.getFolder("INBOX");

        inbox.open(Folder.READ_ONLY);

        return inbox.getMessages();
    }
}