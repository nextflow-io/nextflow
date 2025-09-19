/*
 * Copyright 2024-2025, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

@ScopeName("mail")
@Description("""
    The `mail` scope controls the mail server used to send email notifications.
""")
@CompileStatic
@ToString
class MailConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable Java Mail logging for debugging purposes (default: `false`).
    """)
    final boolean debug

    @ConfigOption
    @Description("""
        Default email sender address.
    """)
    final String from

    @ConfigOption
    @Description("""
    """)
    final Duration sendMailTimeout

    @Description("""
        The `mail.smtp` scope supports any SMTP configuration property in the [Java Mail API](https://javaee.github.io/javamail/).
    
        [Read more](https://javaee.github.io/javamail/docs/api/com/sun/mail/smtp/package-summary.html#properties)
    """)
    final SmtpOpts smtp

    /* required by extension point -- do not remove */
    MailConfig() {}

    MailConfig(Map opts) {
        debug = opts.debug as boolean
        from = opts.from
        sendMailTimeout = opts.sendMailTimeout as Duration
        smtp = opts.smtp != null ? new SmtpOpts(opts.smtp as Map) : null
    }

}

@CompileStatic
class SmtpOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        Host name of the mail server.
    """)
    final String host

    @ConfigOption
    @Description("""
        User password to connect to the mail server.
    """)
    final String password

    @ConfigOption
    @Description("""
        Port number of the mail server.
    """)
    final Integer port

    @ConfigOption
    @Description("""
        User name to connect to the mail server.
    """)
    final String user

    private Map opts

    SmtpOpts(Map opts) {
        host = opts.host
        password = opts.password
        port = opts.port as Integer
        user = opts.user
        this.opts = opts
    }

    Map toMap() {
        return opts
    }

}
