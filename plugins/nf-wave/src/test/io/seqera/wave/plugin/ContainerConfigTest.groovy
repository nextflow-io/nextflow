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
 *
 */

package io.seqera.wave.plugin

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerConfigTest extends Specification {

    @Unroll
    def 'should merge env' () {
        given:
        def config = new ContainerConfig()

        expect:
        config.mergeEnv(ENV1, ENV2) == EXPECTED

        where:
        ENV1                | ENV2              | EXPECTED
        null                | null              | null
        []                  | null              | []
        ['foo=1']           | []                | ['foo=1']
        ['foo=1']           | ['bar=2']         | ['foo=1','bar=2']
        ['foo=1']           | ['foo=2']         | ['foo=2']         // <-- env2 overrides env1
        ['foo=1','baz=3']   | ['foo=2']         | ['foo=2','baz=3'] // <-- overrides 'foo' in env1 and keep as first entry
        ['foo=1']           | ['baz=3','foo=2'] | ['foo=2','baz=3'] // <-- overrides 'foo' in env1 and keep as first entry
    }

    def 'should override entry' () {
        expect:
        new ContainerConfig(entrypoint: LEFT) + new ContainerConfig(entrypoint: RIGHT)
                == new ContainerConfig(entrypoint: EXPECTED)

        where:
        LEFT            | RIGHT         | EXPECTED
        null            | null          | null
        ['entry1.sh']   | null          | ['entry1.sh']
        null            | ['entry2.sh'] | ['entry2.sh']
        ['entry1.sh']   | ['entry2.sh'] | ['entry2.sh']
    }

    def 'should override cmd' () {
        expect:
        new ContainerConfig(cmd: LEFT) + new ContainerConfig(cmd: RIGHT)
                == new ContainerConfig(cmd: EXPECTED)
        where:
        LEFT            | RIGHT         | EXPECTED
        null            | null          | null
        ['cmd1.sh']     | null          | ['cmd1.sh']
        null            | ['cmd2.sh']   | ['cmd2.sh']
        ['cmd1.sh']     | ['cmd2.sh']   | ['cmd2.sh']

    }

    def 'should override workdir' () {
        expect:
        new ContainerConfig(workingDir: LEFT) + new ContainerConfig(workingDir: RIGHT)
                == new ContainerConfig(workingDir: EXPECTED)

        where:
        LEFT            | RIGHT         | EXPECTED
        null            | null          | null
        '/foo'          | null          | '/foo'
        null            | '/bar'        | '/bar'
        '/foo'          | '/bar'        | '/bar'

    }

    def 'should merge env config' () {
        expect:
        new ContainerConfig(env:LEFT) + new ContainerConfig(env: RIGHT) == new ContainerConfig(env: EXPECTED)

        where:
        LEFT                    | RIGHT         | EXPECTED
        null                    | null          | null
        ['alpha=1']             | null          | ['alpha=1']
        null                    |  ['beta=2']   |  ['beta=2']
        ['alpha=1','delta=x']   |   ['beta=2','delta=z']  | ['alpha=1','delta=z','beta=2']
    }

    @Unroll
    def 'should merge layers' () {

        expect:
        new ContainerConfig(layers:LEFT) + new ContainerConfig(layers: RIGHT) == new ContainerConfig(layers: EXPECTED)

        where:
        LEFT                                        | RIGHT     | EXPECTED
        null                                        | null      | null
        []                                          | []        | []
        [new ContainerLayer(location: 'http://x')]  | null      | [new ContainerLayer(location: 'http://x')]
        and:
        [new ContainerLayer(location: 'http://x')]  | [new ContainerLayer(location: 'http://y'),new ContainerLayer(location: 'http://y')]   | [new ContainerLayer(location: 'http://x'), new ContainerLayer(location: 'http://y'),new ContainerLayer(location: 'http://y')]
    }

    def 'should compute config fingerprint' () {
        given:
        def config1 = new ContainerConfig(entrypoint: ['entry.sh'], cmd: ['the.cmd'], env: ['x=2'], workingDir: '/foo')
        def config2 = new ContainerConfig(entrypoint: ['entry.sh'], cmd: ['the.cmd'], env: ['x=2'], workingDir: '/foo')
        def config3 = new ContainerConfig(entrypoint: ['entry.sh'], cmd: ['the.cmd'], env: ['x=2'], workingDir: '/bar')
        and:
        def layer1 = new ContainerLayer('http://this/that', 'abd', 100, 'efg')
        def layer2 = new ContainerLayer('http://this/that', 'abd', 100, 'efg')
        def layer3 = new ContainerLayer('http://xxx/yyy', 'abd', 200, 'efg')
        and:
        def fusion1 = new ContainerLayer('abc', 'xyz', 300, 'efg', true)
        def fusion2 = new ContainerLayer('efg', 'pqr', 400, 'efg', true)

        expect:
        config1.fingerprint() == config2.fingerprint()
        config1.fingerprint() != config3.fingerprint()

        when:
        config1.appendLayer(layer1)
        config2.appendLayer(layer2)
        then:
        // layers have the same fingerprint
        layer1.fingerprint() == layer2.fingerprint()
        // configs have the same fingerprint
        config1.fingerprint() == config2.fingerprint()

        when:
        config1.layers.clear()
        config2.layers.clear()
        and:
        config1.appendLayer(layer1)
        config2.appendLayer(layer3)
        then:
        // different layer fingerprint cause the config to have different fingerprints
        layer1.fingerprint() != layer3.fingerprint()
        config1.fingerprint() != config2.fingerprint()

        when:
        config1.layers.clear()
        config2.layers.clear()
        and:
        config1.appendLayer(fusion1)
        config2.appendLayer(fusion2)
        then:
        // the layer have different fingerprint BUT the `skipHashing`
        // makes the config to ignore it
        fusion1.fingerprint() != fusion2.fingerprint()
        config1.fingerprint() == config2.fingerprint()
    }

    def 'should validate empty' () {
        expect:
        new ContainerConfig().empty()
        new ContainerConfig([], null, null, null, null).empty()
        new ContainerConfig(null, [], null, null, null).empty()
        new ContainerConfig(null, null, [], null, null).empty()
        new ContainerConfig(null, null, null, '', null).empty()
        new ContainerConfig(null, null, null, null, []).empty()
        and:
        !new ContainerConfig(['x'], null, null, null, null).empty()
        !new ContainerConfig(null, ['x'], null, null, null).empty()
        !new ContainerConfig(null, null, ['x'], null, null).empty()
        !new ContainerConfig(null, null, null, 'x', null).empty()
        !new ContainerConfig(null, null, null, null, [new ContainerLayer()]).empty()
    }
}
