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

package nextflow.cli

import spock.lang.Specification
/**
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class CmdPluginTest extends Specification {

    private String capture(List<String> args) {
        final cmd = new CmdPlugin(args: args)
        final buffer = new ByteArrayOutputStream()
        final origOut = System.out
        System.setOut(new PrintStream(buffer))
        try {
            cmd.usage(args)
        }
        finally {
            System.setOut(origOut)
        }
        return buffer.toString()
    }

    def 'should print the general usage when no args are given' () {
        when:
        def out = capture([])
        then:
        out.contains('Usage: nextflow plugin <sub-command> [options]')
        out.contains('install <pluginId,...>')
        out.contains('create [<name> <provider> [dir]]')
        out.contains('<plugin-name>:<command> [args]')
    }

    def 'should print the install usage' () {
        when:
        def out = capture(['install'])
        then:
        out.contains('Usage: nextflow plugin install <pluginId,...>')
    }

    def 'should print the create usage' () {
        when:
        def out = capture(['create'])
        then:
        out.contains('Usage: nextflow plugin create [<plugin name> <provider name> [project path]]')
    }

    def 'should print the plugin-specific command usage for the colon form' () {
        when:
        def out = capture(['nf-hello:greet'])
        then:
        out.contains('Usage: nextflow plugin <plugin-name>:<command> [args]')
        and: 'it does not abort claiming the command is unknown'
        !out.contains('Unknown plugin sub-command')
    }

    def 'should print the failure reason followed by the general usage for an unknown sub-command' () {
        when:
        def out = capture(['instal'])
        then:
        out.contains('Unknown plugin sub-command: instal')
        and: 'the reason is printed before the general usage'
        out.indexOf('Unknown plugin sub-command: instal') < out.indexOf('Usage: nextflow plugin <sub-command> [options]')
    }

}
