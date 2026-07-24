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

package test

import spock.lang.Specification

/**
 * Tests for {@link TestHelper#filterLogNoise}.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class TestHelperTest extends Specification {

    def 'should keep command output and drop log noise'() {
        given:
        def output = '''\
            10:00:00.000 [main] DEBUG nextflow.Session - some debug
            10:00:00.001 [main] INFO  nextflow.Session - some info
            real output line 1
            real output line 2'''.stripIndent()

        when:
        def lines = TestHelper.filterLogNoise(output)

        then:
        lines == ['real output line 1', 'real output line 2']
    }

    def 'should drop an exception stack trace leaking from plugin loading'() {
        given: 'command output polluted by a pf4j plugin-descriptor error stack trace'
        def output = '''\
            No workflow runs found in lineage history log
            org.pf4j.InvalidPluginDescriptorException: Field 'id' cannot be empty
            \tat org.pf4j.AbstractPluginManager.validatePluginDescriptor(AbstractPluginManager.java:990)
            \tat org.pf4j.AbstractPluginManager.loadPluginFromPath(AbstractPluginManager.java:896)
            \tat org.pf4j.AbstractPluginManager.loadPlugins(AbstractPluginManager.java:246)
            Caused by: java.lang.IllegalStateException: boom
            \t... 12 more'''.stripIndent()

        when:
        def lines = TestHelper.filterLogNoise(output)

        then: 'only the real command output remains'
        lines == ['No workflow runs found in lineage history log']
    }

    def 'should not treat ordinary error messages as stack-trace noise'() {
        given:
        def output = 'Error loading lid://12345 - Lineage record 12345 not found'

        when:
        def lines = TestHelper.filterLogNoise(output)

        then:
        lines == ['Error loading lid://12345 - Lineage record 12345 not found']
    }
}
