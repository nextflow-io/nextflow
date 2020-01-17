/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.nio.file.Paths
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId

import nextflow.NextflowMeta
import nextflow.mail.Attachment
import nextflow.mail.Mail
import nextflow.mail.Mailer
import nextflow.trace.WorkflowStats
import nextflow.util.Duration
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowNotifierTest extends Specification {

    def 'should load notification text template' () {

        given:
        def work = TestHelper.createInMemTempDir()
        def sessionId = UUID.randomUUID()
        def now = OffsetDateTime.ofInstant(Instant.ofEpochMilli(1513285947928), ZoneId.systemDefault())
        def end = now.plusSeconds(150)
        def meta = new WorkflowMetadata(
                runName: 'foo_bartali',
                exitStatus: 0,
                start: now,
                complete: end,
                duration: Duration.between( now, end ),
                commandLine: 'nextflow run big-workflow',
                launchDir: Paths.get('/launch/dir'),
                workDir: work,
                projectDir: Paths.get('/project/dir'),
                scriptFile: Paths.get('/script/file.nf'),
                scriptName: 'script.nf',
                scriptId: '9ds9ds090',
                sessionId: sessionId,
                repository: 'http://github.com/big/workflow',
                revision: 'master',
                commitId: '589294',
                profile: 'my-cluster',
                container: 'image/foo:tag',
                containerEngine: 'docker',
                nextflow: new NextflowMeta('0.27.0', 333, '2017-12-12'),
                stats: new WorkflowStats(succeedMillis: 4_000_000, succeededCount: 10, failedCount: 20, cachedCount: 30, ignoredCount: 0)
        )

        def notifier = new WorkflowNotifier(variables: [workflow:meta], config: [:])

        when:
        meta.success = true
        def template = notifier.loadDefaultTextTemplate()
        then:
        template == """
                ====================================
                 Workflow completion notification
                ====================================
                 Run Name: foo_bartali

                ## Execution completed successfully! ##

                The command used to launch the workflow was as follows:

                  nextflow run big-workflow

                ** Execution summary **

                  Launch time       : ${now.format('dd-MMM-yyyy HH:mm:ss')}
                  Ending time       : ${end.format('dd-MMM-yyyy HH:mm:ss')} (duration: 2m 30s)
                  Total CPU-Hours   : 1.1
                  Tasks stats       : Succeeded 10; Cached 30; Ignored 0; Failed 20
                  Launch directory  : /launch/dir
                  Work directory    : ${work.toUriString()}
                  Project directory : /project/dir
                  Script name       : script.nf
                  Script ID         : 9ds9ds090
                  Workflow session  : ${sessionId}
                  Workflow repo     : http://github.com/big/workflow
                  Workflow revision : master (589294)
                  Workflow profile  : my-cluster
                  Workflow container: image/foo:tag
                  Container engine  : docker
                  Nextflow version  : 0.27.0, build 333 (2017-12-12)

                --
                This email was sent by Nextflow
                cite doi:10.1038/nbt.3820
                http://nextflow.io
                """.stripIndent().leftTrim()


        when:
        meta.success = false
        template = notifier.loadDefaultTextTemplate()
        then:
        template == """
                ====================================
                 Workflow completion notification
                ====================================
                 Run Name: foo_bartali

                #########################################
                ## Execution completed unsuccessfully! ##
                #########################################
                The exit status of the task that caused the workflow execution to fail was: 0.
                The full error message was:

                n/a


                The command used to launch the workflow was as follows:

                  nextflow run big-workflow

                ** Execution summary **

                  Launch time       : ${now.format('dd-MMM-yyyy HH:mm:ss')}
                  Ending time       : ${end.format('dd-MMM-yyyy HH:mm:ss')} (duration: 2m 30s)
                  Total CPU-Hours   : 1.1
                  Tasks stats       : Succeeded 10; Cached 30; Ignored 0; Failed 20
                  Launch directory  : /launch/dir
                  Work directory    : ${work.toUriString()}
                  Project directory : /project/dir
                  Script name       : script.nf
                  Script ID         : 9ds9ds090
                  Workflow session  : ${sessionId}
                  Workflow repo     : http://github.com/big/workflow
                  Workflow revision : master (589294)
                  Workflow profile  : my-cluster
                  Workflow container: image/foo:tag
                  Container engine  : docker
                  Nextflow version  : 0.27.0, build 333 (2017-12-12)

                --
                This email was sent by Nextflow
                cite doi:10.1038/nbt.3820
                http://nextflow.io
                """.stripIndent().leftTrim()

    }


    def 'should load notification html template' () {

        given:
        def sessionId = UUID.randomUUID()
        def now = OffsetDateTime.ofInstant(Instant.ofEpochMilli(1513285947928), ZoneId.systemDefault())
        def end = now.plusSeconds(150)
        def workDir = TestHelper.createInMemTempDir()

        def meta = new WorkflowMetadata(
                runName: 'foo_bartali',
                exitStatus: 0,
                start: now,
                complete: end,
                duration: Duration.between( now, end ),
                commandLine: 'nextflow run big-workflow',
                launchDir: Paths.get('/launch/dir'),
                workDir: workDir,
                projectDir: Paths.get('/project/dir'),
                scriptFile: Paths.get('/script/file.nf'),
                scriptName: 'script.nf',
                scriptId: '9ds9ds090',
                sessionId: sessionId,
                repository: 'http://github.com/big/workflow',
                revision: 'master',
                commitId: '589294',
                profile: 'my-cluster',
                container: 'image/foo:tag',
                containerEngine: 'docker',
                nextflow: new NextflowMeta('0.27.0', 333, '2017-12-12'),
                stats: new WorkflowStats(succeedMillis: 4000)
        )

        def notifier = new WorkflowNotifier(workflow: meta, config: [:], variables: [workflow:meta])

        when:
        meta.success = true
        def template = notifier.loadDefaultHtmlTemplate()
        then:
        template.contains('<h1>Workflow completion notification</h1>')
        template.contains('Execution completed successfully')
        template.contains(workDir.toUriString())

        when:
        meta.success = false
        template = notifier.loadDefaultHtmlTemplate()
        then:
        template.contains('<h1>Workflow completion notification</h1>')
        template.contains('Execution completed unsuccessfully')
        template.contains(workDir.toUriString())

    }

    def 'should normalise template list' () {
        given:
        def notifier = new WorkflowNotifier()

        expect:
        notifier.normaliseTemplate0('foo', []) == [new File('foo')]
        notifier.normaliseTemplate0([Paths.get('yyy'), null], [new File('xxx')]) == [new File('xxx'), new File('yyy')]
    }

    def 'should create notification mail' () {

        given:
        Mail mail
        def workflow = new WorkflowMetadata()
        def notifier = Spy(WorkflowNotifier)
        notifier.workflow = workflow
        def attach = Mock(Attachment)

        /*
         * create success completion  *default* notification email
         */
        when:
        workflow.success = true
        workflow.runName = 'foo'
        mail = notifier.createMail([to:'paolo@yo.com', from:'bot@nextflow.com'])
        then:
        1 * notifier.loadDefaultTextTemplate() >> 'TEXT template'
        1 * notifier.loadDefaultHtmlTemplate() >> 'HTML template'
        1 * notifier.loadDefaultLogo() >> attach
        mail.to == 'paolo@yo.com'
        mail.from == 'bot@nextflow.com'
        mail.subject == 'Workflow completion [foo] - SUCCEED'
        mail.text == 'TEXT template'
        mail.body == 'HTML template'
        mail.attachments == [attach]

        /*
         * create fail completion *custom* notification email
         */
        when:
        workflow.success = false
        workflow.runName = 'bar'
        mail = notifier.createMail([to:'alpha@dot.com', from:'beta@dot.com', template: ['/some/file.txt', '/other/file.html'], binding: [one:1, two:2]])
        then:
        1 * notifier.loadMailTemplate(new File('/some/file.txt'), [one:1, two:2]) >> 'TEXT template'
        1 * notifier.loadMailTemplate(new File('/other/file.html'), [one:1, two:2]) >> 'HTML template'
        mail.to == 'alpha@dot.com'
        mail.from == 'beta@dot.com'
        mail.subject == 'Workflow completion [bar] - FAILED'
        mail.body == 'HTML template'
        mail.text == 'TEXT template'

    }


    def 'should send notification email' () {

        given:
        def workflow = Mock(WorkflowMetadata)
        def notifier = Spy(WorkflowNotifier)
        def CONFIG_NOTIFIER = [
            enabled: true,
            from: 'me@nextflow.io',
            to: 'you@nextflow.io'
        ]
        def CONFIG_MAIL = [
                smtp: ['host': 'gmail.com', port: 485]
        ]

        def MAIL = Mock(Mail)
        def MAILER = Mock(Mailer)

        when:
        notifier.workflow = workflow
        notifier.config = [notification: CONFIG_NOTIFIER, mail: CONFIG_MAIL]
        notifier.sendNotification()
        then:
        1 * notifier.createMail(CONFIG_NOTIFIER) >> MAIL
        1 * notifier.createMailer(CONFIG_MAIL) >> MAILER
        1 * MAILER.send(MAIL)

        when: '''
        `enabled` flag is false, notification is NOT sent
        '''
        notifier.config = [notification: [enabled: false, to:'you@dot.com']]
        notifier.sendNotification()
        then:
        0 * notifier.createMail(_) >> null
        0 * notifier.createMailer(_) >> MAILER
        0 * MAILER.send(MAIL)

        when: '''
        notification is implicitly enabled if recipient field is provided
        '''
        notifier.config = [notification: [to:'you@dot.com']]
        notifier.sendNotification()
        then:
        0 * notifier.createMail(_) >> null
        0 * notifier.createMailer(_) >> MAILER
        0 * MAILER.send(MAIL)


    }


}
