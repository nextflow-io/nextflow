/*
 * Copyright 2013-2026, Seqera Labs
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

import static test.ScriptHelper.*

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.read.ListAppender
import nextflow.extension.Bolts
import org.slf4j.LoggerFactory
import test.Dsl2Spec

/**
 * A script declaring one or more agents must warn -- once -- that agents are
 * a preview feature.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentPreviewWarnTest extends Dsl2Spec {

    ListAppender<ILoggingEvent> appender
    Logger logger

    def setup() {
        // the `warn1` dedup cache is JVM-global and throttles for a minute, so it would
        // swallow the warning if another test in the same JVM already emitted it
        clearLoggerCache()
        logger = (Logger) LoggerFactory.getLogger(BaseScript)
        logger.setLevel(Level.WARN)
        appender = new ListAppender<ILoggingEvent>()
        appender.start()
        logger.addAppender(appender)
    }

    def cleanup() {
        logger.detachAppender(appender)
        clearLoggerCache()
    }

    private static void clearLoggerCache() {
        final field = Bolts.getDeclaredField('LOGGER_CACHE')
        field.setAccessible(true)
        (field.get(null) as Map).clear()
    }

    private List<String> warnings() {
        return appender.list*.formattedMessage
    }

    def 'should warn that the agent construct is experimental' () {
        when:
        loadScript(module: true, '''
            nextflow.enable.types = true

            agent alpha {
                model 'openai/gpt-4o'

                input:
                question: String

                output:
                answer: String

                prompt:
                """
                ${question}
                """
            }

            agent beta {
                model 'openai/gpt-4o'

                input:
                question: String

                output:
                answer: String

                prompt:
                """
                ${question}
                """
            }
            ''')

        then: 'the warning is emitted once, not once per agent definition'
        warnings().findAll { it.startsWith('Agents are') } == [
            'Agents are a preview feature -- syntax and behavior may change in future releases'
        ]
    }

    def 'should not warn when the script declares no agent' () {
        when:
        loadScript(module: true, '''
            nextflow.enable.types = true

            process foo {
                output:
                out: String

                script:
                """
                echo hello
                """
            }
            ''')

        then:
        warnings().findAll { it.startsWith('Agents are') } == []
    }

}
