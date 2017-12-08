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

package nextflow.mail
import javax.activation.DataHandler
import javax.activation.FileDataSource
import javax.mail.Message
import javax.mail.MessagingException
import javax.mail.Session
import javax.mail.internet.InternetAddress
import javax.mail.internet.MimeBodyPart
import javax.mail.internet.MimeMessage
import javax.mail.internet.MimeMultipart
import java.nio.file.Path

import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.util.Duration
/**
 * This class implements the send mail functionality
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Edgar Garriga <edgano@@gmail.com>
 */
@Slf4j
class Mailer {

    private long SEND_MAIL_TIMEOUT = 15_000

    private String protocol = 'smtp'

    /**
     * Mail attachments
     */
    private List<File> attachments

    /**
     * Mail session object
     */
    private Session session

    /**
     * Holds mail settings and configuration attributes
     */
    private Map config = [:]

    private Map env = System.getenv()

    Mailer setConfig( Map params ) {
        if( params )
            this.config.putAll(params)

        if( params.attach )
            setAttachment(params.attach)

        return this
    }

    /**
     * Mail TO recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    String getTo() { config.to }

    Mailer setTo(String value) {
        config.to = value
        return this
    }

    /**
     * Mail sender (addresses must follow RFC822 syntax)
     * Multiple addresses can be separated by a comma
     */
    String getFrom() { config.from }

    Mailer setFrom(String value) {
        config.from = value
        return this
    }

    /**
     * Mail subject
     */
    String getSubject() { config.subject }

    Mailer setSubject(String value) {
        config.subject = value
        return this
    }

    /**
     * Mail content
     */
    String getContent() { config.content }

    Mailer setContent(String value) {
        config.content = value
        return this
    }

    /**
     * Mail content mime-type
     */
    String getType() { config.type ?: 'text/plain' }

    Mailer setType(String value) {
        config.type = value
        return this
    }

    /**
     * Mail CC recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    String getCc() { config.cc }

    Mailer setCc(String value) {
        config.cc = value
        return this
    }

    /**
     * Mail CCC recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    String getBcc() { config.bcc }

    Mailer setBcc(String value) {
        config.bcc = value
        return this
    }

    @Memoized
    static protected String findSysMailer() {

        // first try `sendmail`
        if( runCommand("command -v sendmail &>/dev/null") == 0 )
            return 'sendmail'

        else if( runCommand("command -v mail &>/dev/null") == 0  )
            return 'mail'

        log.warn "Missing system mail command -- Make sure to have either `sendmail` or `mail` tool installed otherwise configure your mail server properties in the nextflow config file"
        return null
    }

    static protected int runCommand(String cmd) {
        def proc = ['bash','-c',cmd].execute()
        proc.waitForOrKill(1_000)
        return proc.exitValue()
    }

    /**
     * Get the properties of the system and insert the properties needed to the mailing procedure
     */
    protected Properties createProps() {

        if( config.smtp instanceof Map ) {
            def cfg = [mail: [smtp: config.smtp]] as ConfigObject
            def props = cfg.toProperties()
            props.setProperty('mail.transport.protocol', 'smtp')
            // -- debug for debugging
            log.trace "Mail session properties:\n${dumpProps(props)}"
            return props
        }

        return new Properties()
    }

    private String dumpProps(Properties props) {
        def dump = new StringBuilder()
        props.each {
            if( it.key.toString().contains('password') )
                dump << "  $it.key=xxx\n"
            else
                dump << "  $it.key=$it.value\n"
        }

        dump.toString()
    }

    /**
     * @return The mail {@link Session} object given the current configuration
     */
    protected Session getSession() {
        if( !session ) {
            session = Session.getInstance(createProps())
        }

        return session
    }

    String getHost() {
        getConfig('host')
    }

    protected int getPort() {
        def port = getConfig('port')
        port ? port as int : -1
    }

    protected String getUser() {
        getConfig('user')
    }

    protected String getPassword() {
       getConfig('password')
    }

    protected getConfig( String name ) {
        def key = "${protocol}.${name}"
        def value = config.navigate(key)
        if( !value ) {
            // fallback on env properties
            value = env.get("NXF_${key.toUpperCase().replace('.','_')}".toString())
        }
        return value
    }

    protected void sendViaJavaMail(MimeMessage message) {

        final transport = getSession().getTransport()
        transport.connect(host, port as int, user, password)
        try {
            transport.sendMessage(message, message.getAllRecipients())
        }
        finally {
            transport.close()
        }
    }

    protected long getSendTimeout() {
        def timeout = config.sendMailTimeout as Duration
        return timeout ? timeout.toMillis() : SEND_MAIL_TIMEOUT
    }

    protected void sendViaSendmail(MimeMessage message) {

        def cmd = ['sendmail','-t']
        def proc = new ProcessBuilder()
                        .command(cmd)
                        .redirectErrorStream(true)
                        .start()
        // pipe the message to the sendmail stdin
        final stdout = new StringBuilder()
        final stdin = proc.getOutputStream()
        message.writeTo(stdin);
        stdin.close()   // <-- don't forget otherwise it hangs
        // wait for the sending to complete
        proc.consumeProcessOutputStream(stdout)
        proc.waitForOrKill(sendTimeout)
        def status = proc.exitValue()
        if( status != 0 ) {
            throw new MessagingException("Unable to send mail message\n  sendmail exit status: $status\n  reported error: $stdout")
        }
    }

    /**
     * @return A multipart mime message representing the mail message to send
     */
    protected MimeMessage createMimeMessage() {

        final msg = new MimeMessage(getSession())

        if( subject )
            msg.setSubject(subject)

        if( to )
            msg.setRecipients(Message.RecipientType.TO, to)

        if( cc )
            msg.setRecipients(Message.RecipientType.CC, cc)

        if( bcc )
            msg.setRecipients(Message.RecipientType.BCC, bcc)

        if( from )
            msg.addFrom(InternetAddress.parse(from))

        MimeMultipart multipart = new MimeMultipart()

        if( content ) {
            def main = new MimeBodyPart()
            main.setContent( content, type )
            multipart.addBodyPart(main)
            msg.setContent(multipart)
        }

        // -- attachment
        for( File file : attachments ) {
            def attach = new MimeBodyPart()
            def source = new FileDataSource(file)
            attach.setDataHandler(new DataHandler(source))
            attach.setFileName(file.getName())
            multipart.addBodyPart(attach)
        }

        return msg
    }

    List<File> getAttachments() { attachments }

    Mailer addAttachment( attach ) {
        if( attachments==null )
            attachments = []

        if( attach instanceof Collection ) {
            attach.each { attach0(it) }
        }
        else if( attach instanceof Object[] ) {
            attach.each { attach0(it) }
        }
        else
            attach0(attach)

        return this
    }

    Mailer setAttachment( attach ) {
        attachments = []
        addAttachment(attach)
    }

    /**
     * Add a single item to the list of file attachments
     * @param attach either a {@link File}, {@link Path} or a string file path
     */
    private void attach0( attach ) {
        if( attach instanceof File ) {
            attachments.add(attach)
        }
        else if( attach instanceof Path ) {
            attachments.add(attach.toFile())
        }
        else if( attach instanceof String || attach instanceof GString ) {
            attachments.add(new File(attach.toString()))
        }
        else if( attach != null )
            throw new IllegalArgumentException("Invalid attachment argument: $attach [${attach.getClass()}]")
    }

    /**
     * Send the mail given the provided config setting
     */
    void send() {
        log.trace "Mailer config: $config"

        def msg = createMimeMessage()

        // if the user provided required configuration
        // send via Java Mail API
        if( config.containsKey(protocol) ) {
            log.trace "Mailer send via javamail"
            sendViaJavaMail(msg)
        }
        // otherwise fallback on system sendmail
        else {
            log.trace "Mailer send via sendmail"
            sendViaSendmail(msg)
        }
    }


    /**
     * Send a mail message
     *
     * @param params
     *      The following named parameter can be specified:
     *      - from: mail sender address (multiple addresses can be separated by a comma)
     *      - to: mail TO recipients (multiple addresses can be separated by a comma)
     *      - cc: mail CC recipients (multiple addresses can be separated by a comma)
     *      - bcc: mail BCC recipients (multiple addresses can be separated by a comma)
     *      - content: mail content
     *      - type: mail content mime-type
     *      - attach: file attachment (multiple attachments can be specified as a List of file objects)
     */
    void send( Map params ) {
        config.putAll(params)
        if( params.attach )
            setAttachment(params.attach)

        send()
    }


    /**
     * Send a mail message using a closure to fetch the required parameters
     *
     * @param params
     *    A closure representing the mail message to send eg
     *    <code>
     *        sendMail {
     *          to 'me@dot.com'
     *          from 'your@name.com'
     *          attach '/some/file/path'
     *          subject 'Hello'
     *          content '''
     *           Hi there,
     *           Hope this email find you well
     *          '''
     *        }
     *    <code>
     */
    void send( Closure params ) {
        def wrapper = new MailParams()
        def copy = (Closure)params.clone();
        copy.setResolveStrategy(Closure.DELEGATE_ONLY);
        copy.setDelegate(wrapper);
        copy.call(wrapper);
        send( wrapper.delegate )
    }

}