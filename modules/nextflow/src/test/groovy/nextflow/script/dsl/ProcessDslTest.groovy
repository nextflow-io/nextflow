/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.script.dsl

import spock.lang.Specification

import nextflow.script.BaseScript
import nextflow.script.params.FileInParam
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.TokenVar
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDslTest extends Specification {

    def createDsl() {
        new ProcessDsl(Mock(BaseScript), null)
    }

    def 'should create input directives' () {

        setup:
        def dsl = createDsl()
        def config = dsl.getConfig()

        when:
        dsl._in_file([infile:'filename.fa'])
        dsl._in_val('x')
        dsl._in_stdin()

        then:
        config.getInputs().size() == 3

        config.getInputs().get(0) instanceof FileInParam
        config.getInputs().get(0).name == 'infile'
        (config.getInputs().get(0) as FileInParam).filePattern == 'filename.fa'

        config.getInputs().get(1) instanceof ValueInParam
        config.getInputs().get(1).name == 'x'

        config.getInputs().get(2).name == '-'
        config.getInputs().get(2) instanceof StdInParam

        config.getInputs().names == [ 'infile', 'x', '-' ]
        config.getInputs().ofType( FileInParam ) == [ config.getInputs().get(0) ]

    }

    def 'should create output directives' () {

        setup:
        def dsl = createDsl()
        def config = dsl.getConfig()

        when:
        dsl._out_stdout()
        dsl._out_file(new TokenVar('file1'))
        dsl._out_file(new TokenVar('file2'))
        dsl._out_file(new TokenVar('file3'))

        then:
        config.getOutputs().size() == 4
        config.getOutputs().names == ['-', 'file1', 'file2', 'file3']
        config.getOutputs().ofType(StdOutParam).size() == 1

        config.getOutputs()[0] instanceof StdOutParam
        config.getOutputs()[1].name == 'file1'
        config.getOutputs()[2].name == 'file2'
        config.getOutputs()[3].name == 'file3'

    }

    def 'should clone config object' () {

        when:
        def dsl = createDsl()
        def config = dsl.getConfig()
        dsl.queue 'cn-el6'
        dsl.container 'ubuntu:latest'
        dsl.memory '10 GB'
        dsl._in_val('foo')
        dsl._in_file('sample.txt')
        dsl._out_file('result.txt')

        then:
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().size() == 2
        config.getOutputs().size() == 1

        when:
        def copy = config.clone()
        def builder = new ProcessConfigBuilder(copy)
        builder.queue 'long'
        builder.container 'debian:wheezy'
        builder.memory '5 GB'

        then:
        copy.queue == 'long'
        copy.container == 'debian:wheezy'
        copy.memory == '5 GB'
        copy.getInputs().size() == 2
        copy.getOutputs().size() == 1

        // original config is not affected
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().size() == 2
        config.getOutputs().size() == 1
    }

}
