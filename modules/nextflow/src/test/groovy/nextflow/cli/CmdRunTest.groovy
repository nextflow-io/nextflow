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
 */

package nextflow.cli


import java.nio.file.Files

import nextflow.NextflowMeta
import nextflow.SysEnv
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import spock.lang.Unroll
import test.OutputCapture

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdRunTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @Unroll
    def 'should parse cmd param=#STR' () {

        expect:
        CmdRun.parseParamValue(STR)  == EXPECTED

        where:
        STR         | EXPECTED
        null        | null
        'true'      | true
        'false'     | false
        'foo'       | 'foo'
        '10'        | 10i
        '-10'       | -10i
        '20.00'     | 20i
        '-20.00'    | -20i
        '3000000000'| 3000000000l
        '20.33'     | 20.33d
        '-20.33'    | -20.33d
        '-foo'      | '-foo'
        '--foo'     | '--foo'
        '20x0'      | '20x0'
        '20.d'      | '20.d'
        '20d'       | '20d'
        '20..0'     | '20..0'
        '20..'      | '20..'
        '..20'      | '..20'
    }

    def 'should not detect params type' () {
        given:
        SysEnv.push(NXF_DISABLE_PARAMS_TYPE_DETECTION: 'true')

        expect:
        CmdRun.parseParamValue('true')  == 'true'
        CmdRun.parseParamValue('1000')  == '1000'
        CmdRun.parseParamValue('hola')  == 'hola'

        cleanup:
        SysEnv.pop()
    }

    def 'should parse nested params' () {
        when:
        CmdRun.addParam(PARAMS, KEY, VALUE)
        then:
        PARAMS == EXPECTED

        where:
        PARAMS          | KEY       | VALUE     | EXPECTED
        [:]             | 'foo'     | '1'       | [foo: 1]
        [foo: 1]        | 'bar'     | '2'       | [foo: 1, bar: 2]
        [:]             | 'x.y.z'   | 'Hola'    | [x: [y: [z: 'Hola']]]
        [a: [p:1], x:3] | 'a.q'     | '2'       | [a: [p:1, q: 2], x:3]
        [:]             | /x\.y\.z/ | 'Hola'    | ['x.y.z': 'Hola']
        [:]             | /x.y\.z/  | 'Hola'    | ['x': ['y.z': 'Hola']]
    }

    def 'should convert cli params from kebab case to camel case' () {

        when:
        def params = [:]
        CmdRun.addParam(params, 'alphaBeta', '1')
        CmdRun.addParam(params, 'alpha-beta', '10')
        then:
        params['alphaBeta'] == 10
        !params.containsKey('alpha-beta')

        when:
        params = [:]
        CmdRun.addParam(params, 'aaa-bbb-ccc', '1')
        CmdRun.addParam(params, 'aaaBbbCcc', '10')
        then:
        params['aaaBbbCcc'] == 10
        !params.containsKey('aaa-bbb-ccc')

    }

    def 'should convert kebab case to camel case' () {

        expect:
        CmdRun.kebabToCamelCase('a') == 'a'
        CmdRun.kebabToCamelCase('A') == 'A'
        CmdRun.kebabToCamelCase('a-b-c-') == 'aBC'
        CmdRun.kebabToCamelCase('aa-bb-cc') == 'aaBbCc'
        CmdRun.kebabToCamelCase('alpha-beta-delta') == 'alphaBetaDelta'
        CmdRun.kebabToCamelCase('Alpha-Beta-delta') == 'AlphaBetaDelta'
    }

    @Unroll
    def 'should check run name #STR' () {
        expect:
        CmdRun.matchRunName(STR) == EXPECTED
        where:
        EXPECTED    | STR
        true        | 'foo'
        true        | 'f00'
        true        | 'f-00'
        true        | 'f-a-b'
        true        | 'f-0-1'
        true        | 'foo-bar'
        true        | 'a' * 80
        and:
        true        | 'f_00'
        true        | 'f_a_b'
        true        | 'f_0_1'
        true        | 'foo_bar'
        and:
        false       | '0foo'
        false       | '-foo'
        false       | 'foo--bar'
        false       | 'foo__bar'
        false       | 'foo_-bar'
        false       | 'a' * 81

    }


    def 'should parse params file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def JSON = '{"abc": 1, "xyz": 2}'
        def YAML = '''
                    ---
                    foo: 1
                    bar: 2
                    '''.stripIndent()

        when:
        def file = folder.resolve('params.json')
        file.text = JSON
        and:
        def cmd = new CmdRun(paramsFile: file.toString())
        def params = cmd.parsedParams()
        then:
        params.abc == 1
        params.xyz == 2
        and:
        cmd.hasParams()
        
        when:
        file = folder.resolve('params.yaml')
        file.text = YAML
        and:
        cmd = new CmdRun(paramsFile: file.toString())
        params = cmd.parsedParams()
        then:
        params.foo == 1
        params.bar == 2
        and:
        cmd.hasParams()

        when:
        cmd = new CmdRun(sysEnv: [NXF_PARAMS_FILE: file.toString()])
        params = cmd.parsedParams()
        then:
        params.foo == 1
        params.bar == 2
        and:
        cmd.hasParams()

        when:
        cmd = new CmdRun(sysEnv: [NXF_PARAMS_FILE: '/missing/path.yml'])
        cmd.parsedParams()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Specified params file does not exist: /missing/path.yml'

        cleanup:
        folder?.delete()
    }

    def 'should parse json file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def json = folder.resolve('params.json')
        json.text = '''\
            {
                "alpha": "This is alpha",
                "delta": "${launchDir}/more",
                "gamma": "$should_not_replace",
                "omega": "${baseDir}/end"
            }
            '''.stripIndent()
        when:
        def cmd = new CmdRun(paramsFile: json.toString())
        and:
        def result = cmd.parsedParams( [baseDir: '/BASE/DIR', launchDir: '/WORK/DIR'] )
        then:
        result.alpha == 'This is alpha'
        result.delta == '/WORK/DIR/more'
        result.gamma == '$should_not_replace'
        result.omega == '/BASE/DIR/end'

        cleanup:
        folder.deleteDir()
    }


    def 'should parse yaml file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def json = folder.resolve('params.yaml')
        json.text = '''\
            alpha: "This is alpha"
            delta: 
                beta: "${launchDir}/more"
                gamma: "$should_not_replace"
                omega: "${baseDir}/end"
            '''.stripIndent()
        when:
        def cmd = new CmdRun(paramsFile: json.toString())
        and:
        def result = cmd.parsedParams( [baseDir: '/BASE/DIR', launchDir: '/WORK/DIR'] )
        then:
        result.alpha == 'This is alpha'
        result.delta.beta == '/WORK/DIR/more'
        result.delta.gamma == '$should_not_replace'
        result.delta.omega == '/BASE/DIR/end'

        cleanup:
        folder.deleteDir()
    }

    def 'should check has params' () {
        expect:
        !new CmdRun().hasParams()
        and:
        new CmdRun(params: [foo:'x']).hasParams()
        new CmdRun(paramsFile: '/some/file.yml').hasParams()
        new CmdRun(sysEnv:[NXF_PARAMS_FILE: '/some/file.yml']).hasParams()
    }

    def 'should replace values' () {
        expect:
        // only dollar are ignored
        new CmdRun().replaceVars0('some $xxx there', [xxx:'here']) == 'some $xxx there'

        and:
        // dollar wrapped with {} are interpreted as vars
        new CmdRun().replaceVars0('some ${xxx} ${xxx}', [xxx:'here']) == 'some here here'

        and:
        def text = '''\
            alpha: "${baseDir}/hello"
            delta: "${launchDir}/world"
            gamma: "${012345}"
            omega: "${unknown}"
            '''.stripIndent()
        
        new CmdRun().replaceVars0(text, [baseDir:'/HOME', launchDir: '/WORK' ] ) == '''\
            alpha: "/HOME/hello"
            delta: "/WORK/world"
            gamma: "${012345}"
            omega: "${unknown}"
            '''.stripIndent()
    }

    def 'should validate dont kill jobs' () {
        when:
        def cmd = new CmdRun()
        then:
        cmd.getDisableJobsCancellation() == false

        when:
        cmd = new CmdRun(disableJobsCancellation: true)
        then:
        cmd.getDisableJobsCancellation() == true

        when:
        cmd = new CmdRun(sysEnv: [NXF_DISABLE_JOBS_CANCELLATION: true])
        then:
        cmd.getDisableJobsCancellation() == true
    }

    @Unroll
    def 'should guss is repo' () {
        expect:
        CmdRun.guessIsRepo(PATH) == EXPECTED
        
        where:
        EXPECTED    | PATH
        true        | 'http://github.com/foo'
        true        | 'foo/bar'
        and:
        false       | 'script.nf'
        false       | '/some/path'
        false       | '../some/path'
    }

    def 'should determine dsl mode' () {
        given:
        def DSL1_SCRIPT = '''
        process foo {
          input: 
          file x from ch
        }
        '''

        def DSL2_SCRIPT = '''
        process foo {
          input: 
          file x
        }
        
        workflow { foo() }
        '''

        expect:
        // default to DSL2 if nothing is specified
        CmdRun.detectDslMode(new ConfigMap(), '', [:]) == '2'

        and:
        // take from the config
        CmdRun.detectDslMode(new ConfigMap([nextflow:[enable:[dsl:1]]]), '', [:]) == '1'

        and:
        // the script declaration has priority
        CmdRun.detectDslMode(new ConfigMap([nextflow:[enable:[dsl:1]]]), 'nextflow.enable.dsl=3', [:]) == '3'

        and:
        // env variable is ignored when the config is provided
        CmdRun.detectDslMode(new ConfigMap([nextflow:[enable:[dsl:1]]]), 'echo hello', [NXF_DEFAULT_DSL:'4']) == '1'

        and:
        // env variable is used if nothing else is specified
        CmdRun.detectDslMode(new ConfigMap(), 'echo hello', [NXF_DEFAULT_DSL:'4']) == '4'

        and:
        // dsl mode is taken from the config
        CmdRun.detectDslMode(new ConfigMap([nextflow:[enable:[dsl:4]]]), DSL1_SCRIPT, [:]) == '4'

        and:
        // dsl mode is taken from the config
        CmdRun.detectDslMode(new ConfigMap([nextflow:[enable:[dsl:4]]]), DSL2_SCRIPT, [:]) == '4'

        and:
        // detect version from DSL1 script
        CmdRun.detectDslMode(new ConfigMap(), DSL1_SCRIPT, [NXF_DEFAULT_DSL:'2']) == '1'

        and:
        // detect version from DSL1 script
        CmdRun.detectDslMode(new ConfigMap(), DSL1_SCRIPT, [:]) == '1'

        and:
        // detect version from env
        CmdRun.detectDslMode(new ConfigMap(), DSL2_SCRIPT, [NXF_DEFAULT_DSL:'2']) == '2'

        and:
        // detect version from global default
        CmdRun.detectDslMode(new ConfigMap(), DSL2_SCRIPT, [:]) == '2'
    }

    def 'should warn for invalid config vars' () {
        given:
        def ENV = [NXF_ANSI_SUMMARY: 'true']

        when:
        new CmdRun().checkConfigEnv(new ConfigMap([env:ENV]))

        then:
        def warning = capture
                .toString()
                .readLines()
                .findResults { line -> line.contains('WARN') ? line : null }
                .join('\n')
        and:
        warning.contains('Nextflow variables must be defined in the launching environment - The following variable set in the config file is going to be ignored: \'NXF_ANSI_SUMMARY\'')
    }

    def 'should not warn for valid config vars' () {
        given:
        def ENV = [FOO: '/something', NXF_DEBUG: 'true']

        when:
        new CmdRun().checkConfigEnv(new ConfigMap([env:ENV]))

        then:
        def warning = capture
                .toString()
                .readLines()
                .findResults { line -> line.contains('WARN') ? line : null }
                .join('\n')
        and:
        !warning
    }

    @Unroll
    def 'should detect moduleBinaries' () {
        given:
        NextflowMeta.instance.moduleBinaries(INITIAL)
        CmdRun.detectModuleBinaryFeature(new ConfigMap(CONFIG))

        expect:
        NextflowMeta.instance.isModuleBinariesEnabled() == EXPECTED

        cleanup:
        NextflowMeta.instance.moduleBinaries(false)

        where:
        INITIAL | CONFIG                                          | EXPECTED
        true    | [nextflow: [enable: [ moduleBinaries: true ]]]  | true
        false   | [nextflow: [enable: [ moduleBinaries: true ]]]  | true
        false   | [nextflow: [enable: [ moduleBinaries: false ]]] | false
        true    | [nextflow: [enable: [ moduleBinaries: false ]]] | true
        false   | [:]                                             | false
        true    | [:]                                             | true
    }

    @Unroll
    def 'should detect strict mode' () {
        given:
        NextflowMeta.instance.strictMode(INITIAL)
        CmdRun.detectStrictFeature(new ConfigMap(CONFIG), ENV)

        expect:
        NextflowMeta.instance.isStrictModeEnabled() == EXPECTED

        cleanup:
        NextflowMeta.instance.strictMode(false)

        where:
        INITIAL | CONFIG                                  | ENV                        | EXPECTED
        true    | [nextflow: [enable: [ strict: true ]]]  | [:]                        | true
        false   | [nextflow: [enable: [ strict: true ]]]  | [:]                        | true
        false   | [nextflow: [enable: [ strict: false ]]] | [:]                        | false
        true    | [nextflow: [enable: [ strict: false ]]] | [:]                        | true
        false   | [:]                                     | [:]                        | false
        true    | [:]                                     | [:]                        | true
        true    | [nextflow: [enable: [ strict: true ]]]  | [NXF_ENABLE_STRICT: true ] | true
        false   | [nextflow: [enable: [ strict: true ]]]  | [NXF_ENABLE_STRICT: true ] | true
        false   | [nextflow: [enable: [ strict: false ]]] | [NXF_ENABLE_STRICT: true ] | false
        true    | [nextflow: [enable: [ strict: false ]]] | [NXF_ENABLE_STRICT: true ] | true
        false   | [:]                                     | [NXF_ENABLE_STRICT: true ] | true
        true    | [:]                                     | [NXF_ENABLE_STRICT: true ] | true

    }
}
