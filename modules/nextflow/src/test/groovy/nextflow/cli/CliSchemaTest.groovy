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

import groovy.json.JsonSlurper
import spock.lang.Specification

/**
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class CliSchemaTest extends Specification {

    private static Map parse(String json) {
        return new JsonSlurper().parseText(json) as Map
    }

    def 'should build the root schema with global options and a command index' () {
        given:
        def commands = [ new CmdInfo(), new CmdRun(), new CmdLineage() ] as List<CmdBase>

        when:
        def schema = parse(CliSchema.root(new CliOptions(), commands))

        then:
        schema.name == 'nextflow'
        schema.path == 'nextflow'
        schema.usage.startsWith('nextflow')

        and: 'global options are reported, but meta options are hidden'
        def opts = schema.params.collectMany { it.opts ?: [] }
        '-config' in opts
        !('-h' in opts)
        !('-help-json' in opts)

        and: 'commands are indexed by name with their description'
        schema.subcommands.size() == 3
        schema.subcommands.containsKey('run')
        schema.subcommands.containsKey('info')
        schema.subcommands.run.help

        and: 'aliases are surfaced when present'
        schema.subcommands.lineage.aliases == ['li']

        and: 'hidden options are reported but flagged'
        def debug = schema.params.find { '-debug' in (it.opts ?: []) }
        debug.hidden == true
    }

    def 'should build a command schema with options and arguments' () {
        when:
        def schema = parse(CliSchema.command(new CmdRun()))

        then:
        schema.name == 'run'
        schema.path == 'nextflow run'
        schema.help
        schema.usage.startsWith('nextflow run')

        and: 'a known option is reported with its metadata'
        def resume = schema.params.find { it.name == 'resume' }
        resume.kind == 'option'
        '-resume' in resume.opts

        and: 'a boolean option is flagged'
        def cache = schema.params.find { it.opts == ['-cache'] }
        cache.is_flag == true

        and: 'positional arguments are reported as arguments'
        schema.params.any { it.kind == 'argument' }

        and: 'meta options are excluded'
        !schema.params.collectMany { it.opts ?: [] }.contains('-help-json')
    }

    def 'should expose command aliases' () {
        when:
        def schema = parse(CliSchema.command(new CmdLineage()))

        then:
        schema.name == 'lineage'
        schema.aliases == ['li']
    }

}
