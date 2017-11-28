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
 * This class implements the send mail functionality
 *
 * @author Edgar Garriga <edgano@@gmail.com>
 */
@Slf4j
class SendMail {
    /** Different protocols Nextflow is able to use to send email */
    static final List<String> PROTOCOLS = ['smtp','imap']

    /** the type of the server ie. IMAP or SMTP */
    private String protocol = 'smtp'

    private String host;

    private String port;

    private String user

    private String password

    /**To send the attachment as a file attached or as a text in the body */
    private boolean includeFile

    /** Either system tool `sendmail` or `mail` */
    private String mailer

    private String SYS_MAILER

    {
        SYS_MAILER = findSysMailer()
    }

    SendMail setConfig( Map config ) {
        if( !config )
            return

        if( config.host )
            this.host = config.host
        if( config.port )
            this.port = config.port?.toString()
        if( config.protocol )
            this.protocol = config.protocol
        if( config.user)
            this.user = config.user
        if( config.password )
            this.password = config.password

        return this
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

    /**
     * Get the properties of the system and insert the properties needed to the mailing procedure
     */
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

    /**
     * Nextflow will create a Message to send it using the Java API
     *
     * <li>@param to: The receiver of the email
     * <li>@param from: The sender of the email
     * <li>@param subject: The subject of the email
     * <li>@param body: String with the content of the email
     * <li>@param attachment: Path of the attachment file, can be NULL
     *
     * @return Message Object
     *
     */
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
            // Part two is attachment
            messageBodyPart = new MimeBodyPart();
            messageBodyPart.attachFile(new File(attachment.toString()));
            multipart.addBodyPart(messageBodyPart);
        } else if(attachment && !includeFile){          //send attachment in the body
            messageBodyPart = new MimeBodyPart();
            String fileString = new File(attachment.toString()).text
            messageBodyPart.setText(fileString)
            multipart.addBodyPart(messageBodyPart)
        }
        // Send the complete message parts
        message.setContent(multipart);
        return message
    }

    /**
     * Send the Message using the Java API
     * <li>@param to: The receiver of the email
     * <li>@param from: The sender of the email
     * <li>@param subject: The subject of the email
     * <li>@param body: String with the content of the email
     * <li>@param attachment: Path of the attachment file, can be NULL
     *
     */
    protected void sendByJavaMail(String to, String from, String subject, String body, Path attachment=null) {

        def message = createJavaMessage(to, from, subject, body, attachment)
        // Send message
        Transport.send(message);
    }

    /**
     * Send the Message using the host message service.
     *
     * It could be `sendmail()` or `mail()`
     * To use one or the other, it use @see getMailer()
     *
     * <li>@param to: The receiver of the email
     * <li>@param from: The sender of the email
     * <li>@param subject: The subject of the email
     * <li>@param body: String with the content of the email
     * <li>@param attachment: Path of the attachment file, can be NULL
     *
     * @return Message Object
     *
     */

    // TODO: `from` is not used, looks weird

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

            // TODO: this looks to simple to work, what about the `body`, `from`, `subject` fields?
            // TODO: bash pipe and redirection won't work in this way
            result.add("sendmail")
            result.add("-t")
            result.add(to)
            result.add("<")
            result.add(attachment)
        }
        else if( mailer=='mail' ) {
            if(includeFile){
                //echo "BODY" | mail -s "SUBJECT"  -a /users/cn/egarriga/example.html edgar.garriga@crg.es
                // TODO: as above, echo pipe won't work in this way.
                // The body need to be redirected through the ProcessBuilder. See https://stackoverflow.com/a/28740942/395921
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

    void send( Map params ) {

        final to = params.to as String
        final from = params.from as String
        final subject = params.subject as String
        final body = params.body as String
        final attach = params.attach as Path

        // TODO: it should support `cc` field
        // TODO: `to` and `cc` should be List of email addresses
        send(to, from, subject, body, attach)
    }

}