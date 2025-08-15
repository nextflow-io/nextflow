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

package nextflow.script

import java.nio.file.Path

import groovy.text.GStringTemplateEngine
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.mail.Attachment
import nextflow.util.SysHelper
import nextflow.mail.Mail
import nextflow.mail.Mailer
import nextflow.mail.Notification
/**
 * Send workflow completion notification
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowNotifier {

    /**
     * A map representing the variables defined in the script global scope
     */
    private Map variables

    /**
     * The {@link WorkflowMetadata} object
     */
    private WorkflowMetadata workflow

    WorkflowNotifier(Map variables, WorkflowMetadata workflow) {
        this.variables = variables
        this.workflow = workflow
    }

    /**
     * Send notification email
     *
     * @param config A {@link Map} representing the nextflow configuration object
     */
    void sendNotification(Map config) {
        final notification = new Notification( config.notification as Map ?: Collections.emptyMap() )

        if (!notification.enabled) {
            return
        }

        if (!notification.to) {
            log.warn "Missing notification email target recipients -- Notification is ignored"
            return
        }

        def mail = createMail(notification)
        def mailer = createMailer( config.mail as Map ?: Collections.emptyMap() )
        mailer.send(mail)
    }

    /**
     * Creates {@link Mailer} object that sends the actual email message
     *
     * @param config The {@link Mailer} settings corresponding to the content of the {@code mail} configuration file scope
     * @return A {@link Mailer} object
     */
    protected Mailer createMailer(Map config) {
        return new Mailer(config)
    }

    /**
     * Create notification {@link nextflow.mail.Mail} object given the user parameters
     *
     * @param notification
     */
    protected Mail createMail(Notification notification) {

        def mail = new Mail()
        // -- the subject
        mail.subject("Workflow completion [${workflow.runName}] - ${workflow.success ? 'SUCCEED' : 'FAILED'}")
        // -- the sender
        if (notification.from)
            mail.from(notification.from.toString())
        // -- the recipient(s)
        if (notification.to)
            mail.to(notification.to.toString())

        // -- load user template(s)
        List<File> templates = []
        normaliseTemplate0(notification.template, templates)
        if (templates) {
            final binding = normaliseBindings0(notification.attributes)

            for( File file : templates ) {
                def content = loadMailTemplate(file, binding)
                def plain = file.extension == 'txt'
                if (plain) {
                    mail.text(content)
                } else {
                    mail.body(content)
                }
            }
        }
        // fallback on default embedded template
        else {
            mail.text(loadDefaultTextTemplate())
            mail.body(loadDefaultHtmlTemplate())
            mail.attach(loadDefaultLogo())
        }

        return mail
    }

    protected Map normaliseBindings0(binding) {

        if (binding == null)
            return Map.of()

        if (binding instanceof Map)
            return binding

        throw new IllegalArgumentException("Not a valid binding object [${binding.class.name}]: $binding")
    }

    protected List<File> normaliseTemplate0(entry, List<File> result) {
        if (entry == null)
            return result

        if (entry instanceof File) {
            result << (File) entry
            return result
        }

        if (entry instanceof Path) {
            result << entry.toFile()
            return result
        }

        if (entry instanceof String || entry instanceof GString) {
            result << new File(entry.toString())
            return result
        }

        if (entry instanceof List) {
            entry.each { normaliseTemplate0(it, result) }
            return result
        }

        throw new IllegalArgumentException("Not a valid template file type [${entry.getClass().getName()}]: $entry")
    }

    /**
     * Load notification email template file
     *
     * @return A string representing the email template
     */
    protected String loadMailTemplate(File file, Map binding) {

        if (file && !file.exists())
            throw new IllegalArgumentException("Notification template file does not exist -- check path: $file")

        loadMailTemplate0(file.newInputStream(), binding)
    }

    /**
     * Load and resolve default text email template
     *
     * @return Resolved text template string
     */
    protected String loadDefaultTextTemplate() {
        loadDefaultTemplate0('/nextflow/mail/notification.txt')
    }

    /**
     * Load and resolve default HTML email template
     *
     * @return Resolved HTML template string
     */
    protected String loadDefaultHtmlTemplate() {
        loadDefaultTemplate0('/nextflow/mail/notification.html')
    }

    /**
     * Load the HTML email logo attachment
     * @return A {@link Attachment} object representing the image logo to be included in the HTML email
     */
    protected Attachment loadDefaultLogo() {
        Attachment.resource('/nextflow/mail/nextflow-logo-v2-min.png', contentId: '<nxf-logo>', disposition: 'inline')
    }

    private String loadDefaultTemplate0(String classpathResource) {
        def source = this.class.getResourceAsStream(classpathResource)
        if (!source)
            throw new IllegalArgumentException("Cannot load notification default template -- check classpath resource: $classpathResource")
        loadMailTemplate0(source, [:])
    }

    private String loadMailTemplate0(InputStream source, Map binding) {
        def map = new HashMap()
        if( variables )
            map.putAll(variables)
        if( binding )
            map.putAll(binding)
        // Add SysHelper to template binding so it can be used in templates
        map.put('SysHelper', SysHelper)

        def template = new GStringTemplateEngine().createTemplate(new InputStreamReader(source))
        template.make(map).toString()
    }

}
