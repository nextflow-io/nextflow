/*
 * Copyright 2020-2023, Seqera Labs
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
 *
 */

package nextflow.mail

import javax.mail.MessagingException
import javax.mail.internet.MimeMessage

import groovy.transform.CompileStatic
import nextflow.util.Duration

/**
 * Base class or sending an email via sys command `mail` or `sendmail`
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class BaseMailProvider implements MailProvider {

    private long SEND_MAIL_TIMEOUT = 15_000

    /**
     * Send a email message by using system tool such as `sendmail` or `mail`
     *
     * @param message A {@link MimeMessage} object representing the email to send
     */
    void send(MimeMessage message, Mailer mailer) {
        final cmd = [name(), '-t']
        final proc = new ProcessBuilder()
                .command(cmd)
                .redirectErrorStream(true)
                .start()
        // pipe the message to the sendmail stdin
        final stdout = new StringBuilder()
        final stdin = proc.getOutputStream()
        message.writeTo(stdin);
        stdin.close()   // <-- don't forget otherwise it hangs
        // wait for the sending to complete
        final consumer = proc.consumeProcessOutputStream(stdout)
        proc.waitForOrKill(sendTimeout(mailer))
        def status = proc.exitValue()
        if( status != 0 ) {
            consumer.join()
            throw new MessagingException("Unable to send mail message\n  $mailer exit status: $status\n  reported error: $stdout")
        }
    }

    private long sendTimeout(Mailer mailer) {
        final timeout = mailer.config.sendMailTimeout as Duration
        return timeout ? timeout.toMillis() : SEND_MAIL_TIMEOUT
    }



}
