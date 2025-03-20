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
package nextflow.config.scopes;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class MailConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        When `true` enables Java Mail logging for debugging purpose.
        """)
    public boolean debug;

    @ConfigOption
    @Description("""
        Default email sender address.
        """)
    public String from;

    @Description("""
        The `mail.smtp` scope supports any SMTP configuration property in the [Java Mail API](https://javaee.github.io/javamail/).
    
        [Read more](https://javaee.github.io/javamail/docs/api/com/sun/mail/smtp/package-summary.html#properties)
    """)
    public MailSmtpConfig smtp;

}

class MailSmtpConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Host name of the mail server.
        """)
    public String host;

    @ConfigOption
    @Description("""
        User password to connect to the mail server.
        """)
    public String password;

    @ConfigOption
    @Description("""
        Port number of the mail server.
        """)
    public int port;

    @ConfigOption
    @Description("""
        User name to connect to the mail server.
        """)
    public String user;

}
