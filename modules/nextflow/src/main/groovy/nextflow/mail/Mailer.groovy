/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.mail

import javax.activation.DataHandler
import javax.activation.URLDataSource
import javax.mail.Message
import javax.mail.MessagingException
import javax.mail.Session
import javax.mail.internet.HeaderTokenizer
import javax.mail.internet.InternetAddress
import javax.mail.internet.MimeBodyPart
import javax.mail.internet.MimeMessage
import javax.mail.internet.MimeMultipart
import javax.mail.internet.MimeUtility
import java.nio.charset.Charset
import java.util.regex.Pattern

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.io.LogOutputStream
import nextflow.plugin.Plugins
import org.jsoup.Jsoup
import org.jsoup.nodes.Document
import org.jsoup.parser.Parser
import org.jsoup.safety.Safelist
/**
 * This class implements the send mail functionality
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Edgar Garriga <edgano@gmail.com>
 */
@Slf4j
@CompileStatic
class Mailer {

    // Adapted from post by Phil Haack and modified to match better
    // See https://stackoverflow.com/a/22581832/395921
    private final static String TAG_START = "\\<\\w+((\\s+\\w+(\\s*\\=\\s*(?:\".*?\"|'.*?'|[^'\"\\>\\s]+))?)+\\s*|\\s*)\\>"

    private final static String TAG_END = "\\</\\w+\\>"

    private final static String TAG_SELF_CLOSING = "\\<\\w+((\\s+\\w+(\\s*\\=\\s*(?:\".*?\"|'.*?'|[^'\"\\>\\s]+))?)+\\s*|\\s*)/\\>"

    private final static String HTML_ENTITY = "&[a-zA-Z][a-zA-Z0-9]+;"

    private final static Pattern HTML_PATTERN = Pattern.compile("("+TAG_START+".*"+TAG_END+")|("+TAG_SELF_CLOSING+")|("+HTML_ENTITY+")", Pattern.DOTALL )

    private static String DEF_CHARSET = Charset.defaultCharset().toString()

    /**
     * Mail session object
     */
    private Session session

    /**
     * Holds mail settings and configuration attributes
     */
    private Map config = [:]

    private Map env = System.getenv()

    private String fMailer

    Mailer setConfig( Map params ) {
        if( params ) {
            if( config == null ) config = [:]
            config.putAll(params)
        }
        return this
    }

    Map getConfig() { config }

    protected String getSysMailer() {
        if( !fMailer )
            fMailer = findSysMailer()
        return fMailer
    }

    @Memoized
    static protected String findSysMailer() {

        // first try `sendmail`
        if( runCommand("command -v sendmail &>/dev/null") == 0 )
            return 'sendmail'

        else if( runCommand("command -v mail &>/dev/null") == 0  )
            return 'mail'

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
    @CompileDynamic
    protected Properties createProps() {

        if( config.smtp instanceof Map ) {
            def cfg = [mail: [smtp: config.smtp]] as ConfigObject
            def props = cfg.toProperties()
            props.setProperty('mail.transport.protocol', config.transport?.protocol  ?: 'smtp')
            // -- check proxy configuration
            if( !props.contains('mail.smtp.proxy.host') && System.getProperty('http.proxyHost') ) {
                props['mail.smtp.proxy.host'] = System.getProperty('http.proxyHost')
                props['mail.smtp.proxy.port'] = System.getProperty('http.proxyPort')
            }

            // -- debug for debugging
            if( config.debug == true ) {
                log.debug "Mail session properties:\n${dumpProps(props)}"
            }
            else
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
            if( config.debug != true ) return session

            session.setDebugOut(new PrintStream( new LogOutputStream() {
                @Override protected void processLine(String line, int logLevel) {
                    log.debug(line)
                }
            } ))
            session.setDebug(true)
        }

        return session
    }

    /**
     * @return The SMTP host name or IP address
     */
    protected String getHost() {
        getConfig('host')
    }

    /**
     * @return The SMTP host port
     */
    protected int getPort() {
        def port = getConfig('port')
        port ? port as int : -1
    }

    /**
     * @return The SMTP user name
     */
    protected String getUser() {
        getConfig('user')
    }

    /**
     * @return The SMTP user password
     */
    protected String getPassword() {
        getConfig('password')
    }

    protected getConfig(String name) {
        def key = "smtp.${name}"
        def value = config.navigate(key)
        if( !value ) {
            // fallback on env properties
            value = env.get("NXF_${key.toUpperCase().replace('.','_')}".toString())
        }
        return value
    }

    /**
     * @return A multipart mime message representing the mail message to send
     */
    protected MimeMessage createMimeMessage0(Mail mail) {

        final msg = new MimeMessage(getSession())

        if( mail.subject )
            msg.setSubject(mail.subject)

        if( mail.from )
            msg.addFrom(InternetAddress.parse(mail.from))
        else if( config.from )
            msg.addFrom(InternetAddress.parse(config.from.toString()))

        if( mail.to )
            msg.setRecipients(Message.RecipientType.TO, mail.to)

        if( mail.cc )
            msg.setRecipients(Message.RecipientType.CC, mail.cc)

        if( mail.bcc )
            msg.setRecipients(Message.RecipientType.BCC, mail.bcc)

        return msg
    }


    /**
     * @return A multipart mime message representing the mail message to send
     */
    protected MimeMessage createMimeMessage(Mail mail) {

        final message = createMimeMessage0(mail)

        final wrap = new MimeBodyPart();
        final cover = new MimeMultipart("alternative")
        final charset = mail.charset ?: DEF_CHARSET

        if( mail.text ) {
            def part = new MimeBodyPart()
            part.setText(mail.text, charset)
            cover.addBodyPart(part)
        }

        if( mail.body ) {
            def part = new MimeBodyPart()
            String type = mail.type ?: guessMimeType(mail.body)
            if( !type.contains('charset=') )
                type = "$type; charset=${MimeUtility.quote(charset, HeaderTokenizer.MIME)}"
            part.setContent(mail.body, type)
            cover.addBodyPart(part)
        }
        wrap.setContent(cover)

        // use a separate multipart for body + attachment
        final content = new MimeMultipart("related");
        content.addBodyPart(wrap);

        // -- attachment
        def allFiles = mail.attachments ?: Collections.<Attachment>emptyList()
        for( Attachment item : allFiles ) {
            content.addBodyPart(createAttachment(item))
        }

        message.setContent(content);
        return message
    }

    protected MimeBodyPart createAttachment(Attachment item) {
        final result = new MimeBodyPart()
        if( item.file ) {
            if( !item.file.exists() )
                throw new MessagingException("The following attachment file does not exist: $item.file")
            result.attachFile(item.file)
        }
        else if( item.resource ) {
            def url = this.class.getResource(item.resource)
            if( !url )
                throw new MessagingException("The following attachment resource does not exist: $item.resource")
            def source = new URLDataSource(url)
            result.setDataHandler(new DataHandler(source))
        }
        else {
            throw new IllegalStateException("Invalid attachment object")
        }

        if( item.disposition )
            result.setDisposition(item.disposition)

        if( item.contentId )
            result.setContentID(item.contentId)

        if( item.fileName )
            result.setFileName(item.fileName)

        if( item.description )
            result.setDescription(item.description)

        return result
    }

    /**
     * Creates a pure text email message. It cannot contains attachments
     *
     * @param mail The {@link Mail} object representing the message to send
     * @return A {@link MimeMessage} object instance
     */
    protected MimeMessage createTextMessage(Mail mail) {
        final result = createMimeMessage0(mail)
        final charset = mail.charset ?: DEF_CHARSET
        final text = mail.text ?: stripHtml(mail.body)
        result.setText(text, charset)
        return result
    }

    /**
     * Converts an HTML text to a plain text message
     *
     * @param html The html string to strip
     */
    protected String stripHtml(String html) {
        if( !html )
            return html

        if( !guessHtml(html) )
            return html

        Document document = Jsoup.parse(html);
        document.outputSettings(new Document.OutputSettings().prettyPrint(false));//makes html() preserve linebreaks and spacing
        document.select("br").append("\\n");
        document.select("p").prepend("\\n");
        String s = document.html().replaceAll("\\\\n", "\n");
        def result = Jsoup.clean(s, "", Safelist.none(), new Document.OutputSettings().prettyPrint(false));
        Parser.unescapeEntities(result, false)
    }

    protected boolean guessHtml(String str) {
        str ? HTML_PATTERN.matcher(str).find() : false
    }

    protected String guessMimeType(String str) {
        guessHtml(str) ? 'text/html' : 'text/plain'
    }

    protected boolean detectAwsEnv() {
        if( env.get('AWS_REGION') ) return true
        if( env.get('AWS_DEFAULT_REGION') ) return true
        return false
    }

    protected MailProvider provider() {
        // load all providers
        final providers = Plugins.getExtensions(MailProvider)
        // find the AWS provider
        final awsProvider = providers.find(it -> it.name()=='aws-ses')
        // check if it can use the aws provider
        if( env.get('NXF_ENABLE_AWS_SES')=='true' ) {
            if( awsProvider )
                return awsProvider
            else
                log.warn "Unable to load AWS Simple Email Service (SES) client"
        }

        if( config.containsKey('smtp') ) {
            return providers.find(it -> it.name()=='javamail')
        }

        if( awsProvider && detectAwsEnv() ) {
            return awsProvider
        }

        // detect the mailer type
        final type = getSysMailer()
        return providers.find(it -> it.name()==type)
    }

    /**
     * Send the mail given the provided config setting
     */
    void send(Mail mail) {
        log.trace "Mailer config: $config -- mail: $mail"

        final p = provider()
        if( p != null ) {
            log.debug "Sending mail via `${p.name()}`"
            final msg = p.textOnly()
                    ? createTextMessage(mail)
                    : createMimeMessage(mail)
            p.send(msg, this)
            return
        }

        final msg = "Cannot send email message -- Make sure you have installed `sendmail` or `mail` program or configure a mail SMTP server in the nextflow config file"
        throw new IllegalArgumentException(msg)
    }

    /**
     * Send a mail given a parameter map
     *
     * @param params
     *      The following named parameters are supported
     *      - from: the email sender address
     *      - to: the email recipient address
     *      - cc: the CC recipient address
     *      - bcc: the BCC recipient address
     *      - subject: the email subject
     *      - charset: the email content charset
     *      - type: the email body mime-type
     *      - text: the email plain text alternative content
     *      - body: the email body content (HTML)
     *      - attach: he email attachment
     */
    void send(Map params) {
        send(Mail.of(params))
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
     *          body '''
     *           Hi there,
     *           Hope this email find you well
     *          '''
     *        }
     *    <code>
     */
    void send( Closure params ) {
        def mail = new Mail()
        def copy = (Closure)params.clone()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.setDelegate(mail)
        def body = copy.call(mail)
        if( !mail.body && (body instanceof String || body instanceof GString))
            mail.body(body)
        send(mail)
    }

}
