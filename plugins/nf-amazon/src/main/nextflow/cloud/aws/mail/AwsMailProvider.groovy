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

package nextflow.cloud.aws.mail

import javax.mail.internet.MimeMessage
import software.amazon.awssdk.core.SdkBytes
import software.amazon.awssdk.services.ses.SesClient
import software.amazon.awssdk.services.ses.model.RawMessage
import software.amazon.awssdk.services.ses.model.SendRawEmailRequest
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.mail.MailProvider
import nextflow.mail.Mailer

/**
 * Send a mime message via AWS SES raw API
 *
 * https://docs.aws.amazon.com/ses/latest/dg/send-email-raw.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class AwsMailProvider implements MailProvider {

    @Override
    String name() {
        return 'aws-ses'
    }

    @Override
    boolean textOnly() {
        return false
    }

    @Override
    void send(MimeMessage message, Mailer mailer) {
        final client = getEmailClient()
        // dump the message to a buffer
        final outputStream = new ByteArrayOutputStream()
        message.writeTo(outputStream)
        // send the email
        final rawMessage = RawMessage.builder().data(SdkBytes.fromByteArray(outputStream.toByteArray())).build()
        final result = client.sendRawEmail(SendRawEmailRequest.builder().rawMessage(rawMessage).build())
        log.debug "Mail message sent: ${result}"
    }

    SesClient getEmailClient() {
        return SesClient.builder().build()
    }

}
