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

package nextflow.plugin

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultPluginsTest extends Specification {

    def 'should validate get plugins' () {
        given:
        def AMAZON = new PluginRef('nf-amazon', '0.1.0')
        def GOOGLE = new PluginRef('nf-google', '0.2.0')
        def AZURE = new PluginRef('nf-azure', '0.3.0')
        and:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': AMAZON,
                'nf-google': GOOGLE,
                'nf-azure': AZURE  ])

        expect:
        defaults.getPlugin('nf-amazon') == AMAZON
        defaults.hasPlugin('nf-amazon')
        and:
        !defaults.hasPlugin('nf-foo')

        when:
        defaults.getPlugin('nf-foo')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Unknown Nextflow plugin 'nf-foo'"
    }

    def 'should get sorted string' () {
        given:
        def DELTA = new PluginRef('nf-delta', '1.0.0')
        def OMEGA = new PluginRef('nf-omega', '2.0.0')
        def ALPHA = new PluginRef('nf-alpha', '3.0.0')
        and:
        def defaults = new DefaultPlugins(plugins: [
                'nf-delta': DELTA,
                'nf-omega': OMEGA,
                'nf-alpha': ALPHA  ])

        expect:
        defaults.toSortedString() == 'nf-alpha@3.0.0,nf-delta@1.0.0,nf-omega@2.0.0'
    }

}
