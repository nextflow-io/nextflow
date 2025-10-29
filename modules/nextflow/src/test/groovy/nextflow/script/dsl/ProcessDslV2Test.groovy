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

import java.nio.file.Path

import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessDslV2Test extends Specification {

    def createDsl() {
        new ProcessDslV2(Mock(BaseScript), null)
    }

    def 'should declare process inputs' () {

        given:
        def dsl = createDsl()
        def config = dsl.getConfig()

        when:
        dsl._input_('infile', Path, false)
        dsl._input_('x', String, false)
        dsl._input_('y', String, false)
        dsl.stageAs('filename.fa', { infile })
        dsl.stdin { y }

        then:
        config.getInputs().size() == 3

        config.getInputs().getParams().get(0).name == 'infile'
        config.getInputs().getParams().get(0).type == Path
        config.getInputs().getFiles().get(0).filePattern == 'filename.fa'

        config.getInputs().getParams().get(1).name == 'x'
        config.getInputs().getParams().get(1).type == String

        config.getInputs().getParams().get(2).name == 'y'
        config.getInputs().getParams().get(2).type == String

    }

    def 'should declare process outputs' () {

        given:
        def dsl = createDsl()
        def config = dsl.getConfig()

        when:
        dsl._output_('x', String, { stdout() })
        dsl._output_('file1', Path, { file('$path0') })
        dsl._output_('file2', Path, { file('$path1') })
        dsl._output_('file3', Path, { file('$path2') })
        dsl._unstage_files('$path0', { file1 })
        dsl._unstage_files('$path1', { file2 })
        dsl._unstage_files('$path2', { file3 })

        then:
        config.getOutputs().getParams().size() == 4

        config.getOutputs().getParams().get(0).name == 'x'
        config.getOutputs().getParams().get(1).name == 'file1'
        config.getOutputs().getParams().get(2).name == 'file2'
        config.getOutputs().getParams().get(3).name == 'file3'
        config.getOutputs().getFiles().keySet() == [ '$path0', '$path1', '$path2' ].toSet()
    }

    def 'should declare process topic emissions' () {

        given:
        def dsl = createDsl()
        def config = dsl.getConfig()

        when:
        dsl._topic_({ stdout() }, 'versions')
        dsl._topic_({ file('$path0') }, 'versions')
        dsl._unstage_files('$path0', { file1 })

        then:
        config.getOutputs().getTopics().size() == 2

        config.getOutputs().getTopics().get(0).target == 'versions'
        config.getOutputs().getTopics().get(1).target == 'versions'
        config.getOutputs().getFiles().keySet() == [ '$path0' ].toSet()
    }

    def 'should clone config object' () {

        when:
        def dsl = createDsl()
        def config = dsl.getConfig()
        dsl.queue 'cn-el6'
        dsl.container 'ubuntu:latest'
        dsl.memory '10 GB'
        dsl._input_('foo', String, false)
        dsl._input_('sample', Path, false)
        dsl.stageAs('sample.txt', { sample })
        dsl._output_('result', Path, { file('$file0') })
        dsl._unstage_files('$file0', 'result.txt')

        then:
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().getParams().size() == 2
        config.getOutputs().getParams().size() == 1

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
        copy.getInputs().getParams().size() == 2
        copy.getOutputs().getParams().size() == 1

        // original config is not affected
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().getParams().size() == 2
        config.getOutputs().getParams().size() == 1
    }

}
