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

package nextflow.processor
import java.nio.file.Files
import java.nio.file.Paths

import ch.grengine.Grengine
import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.executor.Executor
import nextflow.file.FileHolder
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.ScriptBinding
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TaskBody
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class TaskRunTest extends Specification {

    def testGetInputsByType() {

        setup:
        def binding = new Binding('x': 1, 'y': 2)
        def task = new TaskRun()
        def list = []

        task.setInput( new StdInParam(binding,list) )
        task.setInput( new FileInParam(binding, list).bind(new TokenVar('x')), 'file1' )
        task.setInput( new FileInParam(binding, list).bind(new TokenVar('y')), 'file2' )
        task.setInput( new EnvInParam(binding, list).bind('z'), 'env' )


        when:
        def files = task.getInputsByType(FileInParam)
        then:
        files.size() == 2

        files.keySet()[0] instanceof FileInParam
        files.keySet()[1] instanceof FileInParam

        files.keySet()[0].name == 'x'
        files.keySet()[1].name == 'y'

        files.values()[0] == 'file1'
        files.values()[1] == 'file2'

    }

    def testGetOutputsByType() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        task.setOutput( new FileOutParam(binding, list).bind('x'), 'file1' )
        task.setOutput( new FileOutParam(binding, list).bind('y'), 'file2' )
        task.setOutput( new StdOutParam(binding, list), 'Hello' )


        when:
        def files = task.getOutputsByType(FileOutParam)
        then:
        files.size() == 2

        files.keySet()[0] instanceof FileOutParam
        files.keySet()[1] instanceof FileOutParam

        files.keySet()[0].filePattern == 'x'
        files.keySet()[1].filePattern == 'y'

        files.values()[0] == 'file1'
        files.values()[1] == 'file2'

    }

    def testGetInputFiles() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        def x = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def y = new FileInParam(binding, list).bind('y')

        task.setInput(x, 1)
        task.setInput(y, [ new FileHolder(Paths.get('file_y_1')) ])

        expect:
        task.getInputFiles().size() == 1
        task.stagedInputs.size() == 1

    }

    def testGetInputFilesMap() {
        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        def x = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def y = new FileInParam(binding, list).bind('y')
        def z = new FileInParam(binding, list).bind('z')

        task.setInput(x, 1)
        task.setInput(y, [ new FileHolder(Paths.get('file_y_1')).withName('foo.txt') ])
        task.setInput(z, [ new FileHolder(Paths.get('file_y_2')).withName('bar.txt') ])

        expect:
        task.getInputFilesMap() == ['foo.txt': Paths.get('file_y_1'), 'bar.txt': Paths.get('file_y_2')]

    }


    def testGetOutputFilesNames() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        when:
        def i1 = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def o1 = new FileOutParam(binding,list).bind('file_out.alpha')
        def o2 = new ValueOutParam(binding,list).bind( 'x' )
        def o3 = new FileOutParam(binding,list).bind('file_out.beta')

        task.setInput(i1, 'Hello' )
        task.setOutput(o1)
        task.setOutput(o2)
        task.setOutput(o3)

        then:
        task.getOutputFilesNames() == ['file_out.alpha', 'file_out.beta']
    }

    def testHasCacheableValues() {

        given:
        def binding = new Binding()
        def list = []

        /*
         * just file as output => no cacheable values
         */
        when:
        def o1 = new FileOutParam(binding,list).bind('file_out.beta')
        def task1 = new TaskRun()
        task1.setOutput(o1)
        then:
        !task1.hasCacheableValues()

        /*
         * when a 'val' is declared as output ==> true
         */
        when:
        def o2 = new ValueOutParam(binding,list).bind( 'x' )
        def task2 = new TaskRun()
        task2.setOutput(o2)
        then:
        task2.hasCacheableValues()

        /*
         * file with parametric name => true
         */
        when:
        def s3 = new FileOutParam(binding, list).bind( new TokenVar('y') )
        def task3 = new TaskRun()
        task3.setOutput(s3)
        then:
        task3.hasCacheableValues()

        when:
        def task5 = new TaskRun( config: new TaskConfig(alpha: 1, beta: 2) )
        then:
        !task5.hasCacheableValues()

        when:
        def task6 = new TaskRun( config: new TaskConfig(alpha: { 'dynamic val' }) )
        then:
        task6.hasCacheableValues()

    }



    def testDumpStdout() {

        setup:
        def file = Files.createTempFile('dumptest',null)
        def task = new TaskRun()

        when:
        task.stdout = '1\n2'
        then:
        task.dumpStdout(5) == ['1','2']

        when:
        task.stdout = """
            1
            2
            3
            4
            5
            6
            7
            8
            9
            """.stripIndent()
        then:
        task.dumpStdout(5) == ['(more omitted..)', '5','6','7','8','9']


        when:
        task = new TaskRun()
        task.stdout = file
        file.text = """
            a
            b
            c
            d
            e
            f
            g
            h
            i
            """.stripIndent()
        then:
        task.dumpStdout(5) == ['e','f','g','h','i']


        when:
        task = new TaskRun()
        task.stdout = Paths.get('/no/file')
        then:
        task.dumpStdout(5) == []

        cleanup:
        file?.delete()

    }

    def testIsSuccess() {

        when:
        def task = new TaskRun(config: [validExitStatus: [0]])
        then:
        task.isSuccess(0) == true
        task.isSuccess('0') == true

        task.isSuccess(1) == false
        task.isSuccess(null) == false
        task.isSuccess('1') == false
        task.isSuccess(Integer.MAX_VALUE) == false


        when:
        task = new TaskRun(config: [validExitStatus: 0])
        then:
        task.isSuccess(0) == true
        task.isSuccess('0') == true

        task.isSuccess(1) == false
        task.isSuccess(null) == false
        task.isSuccess('1') == false
        task.isSuccess(Integer.MAX_VALUE) == false

        when:
        task = new TaskRun(config: [validExitStatus: [0,1]])
        then:
        task.isSuccess(0) == true
        task.isSuccess(1) == true
        task.isSuccess(2) == false

        when:
        task = new TaskRun(config: [validExitStatus: [0]])
        task.exitStatus = 0
        then:
        task.isSuccess() == true

        when:
        task = new TaskRun(config: [validExitStatus: [0]])
        task.exitStatus = 1
        then:
        task.isSuccess() == false

        when:
        task = new TaskRun(config: [validExitStatus: [0]])
        then:
        task.isSuccess() == false
    }

    def 'should return the specified container name' () {
        given:
        def task = [:] as TaskRun
        task.processor = Mock(TaskProcessor)

        when:
        task.config = [container:'foo/bar']
        task.processor.getSession() >> sess
        then:
        task.getContainer() == expected

        where:
        sess                                            | expected
        new Session()                                   | 'foo/bar'
        new Session( docker: [registry:'my.registry'] ) | 'my.registry/foo/bar'

    }

    @Unroll
    def 'should return engine type' () {
        given:
        def task = [:] as TaskRun
        task.processor = Mock(TaskProcessor)
        task.config = new TaskConfig( [container: 'busybox'] )
        task.processor.getSession() >> new Session([(engine): config])

        expect:
        task.container == contnr
        task.containerConfig == config as ContainerConfig
        task.containerConfig.enabled
        task.containerConfig.engine == engine

        where:
        engine      | contnr                    | config
        'docker'    | 'busybox'                 | [enabled: true, x:'alpha', y: 'beta']
        'docker'    | 'd.reg/busybox'           | [enabled: true, x:'alpha', y: 'beta', registry: 'd.reg']
        'udocker'   | 'busybox:latest'          | [enabled: true, x:'alpha', y: 'beta']
        'shifter'   | 'docker:busybox:latest'   | [enabled: true, x:'delta', y: 'gamma']
    }

    def 'should return container image name' () {

        given:
        def task = Spy(TaskRun)
        task.processor = Mock(TaskProcessor)

        when:
        task.script = 'bwa-mem --this'
        task.config = new TaskConfig([container: CONTAINER])
        def image = task.getContainer()
        then:
        task.getContainerConfig() >> [docker:[enabled: true]]
        image == EXPECTED

        where:
        CONTAINER         | EXPECTED
        null              | null
        false             | null
        'debian:latest'   | 'debian:latest'

    }



    def 'should render template and set task attributes'() {

        given:
        // global binding object attached to the script
        def binding = new ScriptBinding([:])
        binding.setParams(query: '/some/file')
        def script = Mock(Script)
        script.getBinding() >> binding

        // local values
        def local = [name: 'Foo']

        // the script template
        def tpl = TestHelper.createInMemTempFile('template.sh')
        tpl.text = 'echo: ${name} ~ query: ${params.query}'

        // the task run
        def task = new TaskRun()
        task.context = new TaskContext(script, local, 'hello')
        task.processor = [:] as TaskProcessor
        task.processor.grengine = new Grengine()

        when:
        def template = task.renderTemplate(tpl)
        then:
        template == 'echo: Foo ~ query: /some/file'
        task.source == 'echo: ${name} ~ query: ${params.query}'
        task.template == tpl

    }

    def 'should resolve the task body script' () {

        given:
        def task = new TaskRun()
        task.processor = [:] as TaskProcessor
        task.processor.grengine = new Grengine()

        /*
         * plain task script
         */
        when:
        task.resolve(new TaskBody({-> 'Hello'}, 'Hello', 'script'))
        then:
        task.script == 'Hello'
        task.source == 'Hello'

        /*
         * task script with a groovy variable
         */
        when:
        task.context = new TaskContext(Mock(Script),[x: 'world'],'foo')
        task.resolve(new TaskBody({-> "Hello ${x}"}, 'Hello ${x}', 'script'))
        then:
        task.script == 'Hello world'
        task.source == 'Hello ${x}'

    }

    def 'should parse a `shell` script' () {

        given:
        // global binding object attached to the script
        def binding = new ScriptBinding([:])
        binding.setParams(var_no: 'NO')
        def script = Mock(Script)
        script.getBinding() >> binding

        def local = [nxf_var: 'YES']

        def task = new TaskRun()
        task.processor = [:] as TaskProcessor
        task.processor.grengine = new Grengine()

        /*
         * bash script where $ vars are ignore and !{xxx} variables are interpolated
        */
        when:
        task.context = new TaskContext(script,local,'foo')
        task.config = new TaskConfig().setContext(task.context)
        task.resolve(new TaskBody({-> '$BASH_VAR !{nxf_var} - !{params.var_no}'}, '<the source script>', 'shell'))  // <-- note: 'shell' type
        then:
        task.script == '$BASH_VAR YES - NO'
        task.source == '<the source script>'

        /*
         * bash script where $ vars are ignore and #{xxx} variables are interpolated
         */
        when:
        task.context = new TaskContext(Mock(Script),[nxf_var: '>interpolated value<'],'foo')
        task.config = new TaskConfig().setContext(task.context)
        task.config.placeholder = '#' as char
        task.resolve(new TaskBody({-> '$BASH_VAR #{nxf_var}'}, '$BASH_VAR #{nxf_var}', 'shell'))  // <-- note: 'shell' type
        then:
        task.script == '$BASH_VAR >interpolated value<'
        task.source == '$BASH_VAR #{nxf_var}'

    }

    def 'should resolve a task template file' () {

        given:
        def task = new TaskRun()
        task.processor = [:] as TaskProcessor
        task.processor.grengine = new Grengine()

        // create a file template
        def file = TestHelper.createInMemTempFile('template.sh')
        file.text = 'echo ${say_hello}'
        // create the task context with two variables
        // - my_file
        // - say_hello
        task.context = new TaskContext(Mock(Script),[say_hello: 'Ciao mondo', my_file: file],'foo')
        task.config = new TaskConfig().setContext(task.context)

        when:
        task.resolve( new TaskBody({-> template(my_file)}, 'template($file)', 'script'))
        then:
        task.script == 'echo Ciao mondo'
        task.source == 'echo ${say_hello}'
        task.template == file

    }

    def 'should resolve a shell template file, ignore BASH variables and parse !{xxx} ones' () {

        given:
        def task = new TaskRun()
        task.processor = [:] as TaskProcessor
        task.processor.grengine = new Grengine()

        // create a file template
        def file = TestHelper.createInMemTempFile('template.sh')
        file.text = 'echo $HOME ~ !{user_name}'
        // create the task context with two variables
        // - my_file
        // - say_hello
        task.context = new TaskContext(Mock(Script),[user_name: 'Foo bar', my_file: file],'foo')
        task.config = new TaskConfig().setContext(task.context)

        when:
        task.resolve( new TaskBody({-> template(my_file)}, 'template($file)', 'shell'))
        then:
        task.script == 'echo $HOME ~ Foo bar'
        task.source == 'echo $HOME ~ !{user_name}'
        task.template == file

    }

    def 'should resolve a shell template file, ignore BASH variables and parse #{xxx} ones' () {

        given:
        def task = new TaskRun()
        task.processor = [:] as TaskProcessor
        task.processor.grengine = new Grengine()

        // create a file template
        def file = TestHelper.createInMemTempFile('template.sh')
        file.text = 'echo $HOME ~ #{user_name}'
        // create the task context with two variables
        // - my_file
        // - say_hello
        task.context = new TaskContext(Mock(Script),[user_name: 'Foo bar', my_file: file],'foo')
        task.config = new TaskConfig().setContext(task.context)
        task.config.placeholder = '#' as char

        when:
        task.resolve( new TaskBody({-> template(my_file)}, 'template($file)', 'shell'))
        then:
        task.script == 'echo $HOME ~ Foo bar'
        task.source == 'echo $HOME ~ #{user_name}'
        task.template == file

    }

    def 'should check container native flag' () {

        given:
        def executor = Mock(Executor)
        def task = Spy(TaskRun);
        task.processor = Mock(TaskProcessor)

        when:
        def result = task.isContainerNative()
        then:
        1 * task.processor.getExecutor() >> executor
        result == false

        when:
        result = task.isContainerNative()
        then:
        1 * task.processor.getExecutor() >> executor
        1 * executor.isContainerNative() >> true
        result == true
    }

    def 'should check container enabled flag' () {

        given:
        def task = Spy(TaskRun);

        when:
        def enabled = task.isContainerEnabled()
        then:
        1 * task.getConfig() >> new TaskConfig()
        !enabled

        when:
        enabled = task.isContainerEnabled()
        then:
        // NO container image is specified => NOT enable even if `enabled` flag is set to true
        _ * task.getConfig() >> new TaskConfig()
        _ * task.getContainerConfig() >> new ContainerConfig([enabled: true])
        _ * task.isContainerNative() >> false
        !enabled

        when:
        enabled = task.isContainerEnabled()
        then:
        // container is specified, not enabled
        _ * task.getConfig() >> new TaskConfig(container:'foo/bar')
        _ * task.getContainerConfig() >> new ContainerConfig([:])
        _ * task.isContainerNative() >> false
        !enabled

        when:
        enabled = task.isContainerEnabled()
        then:
        // container is specified AND native executor (eg kubernetes) => enabled
        _ * task.getConfig() >> new TaskConfig(container:'foo/bar')
        _ * task.getContainerConfig() >> new ContainerConfig([:])
        _ * task.isContainerNative() >> true
        enabled

        when:
        enabled = task.isContainerEnabled()
        then:
        // container is specified AND enabled => enabled
        _ * task.getConfig() >> new TaskConfig(container:'foo/bar')
        _ * task.getContainerConfig() >> new ContainerConfig([enabled: true])
        _ * task.isContainerNative() >> false
        enabled

    }

    def 'should get task environment' () {

        given:
        def EXPECT = [FOO: 'hola', BAR: 'mundo', OMEGA: 'ooo',_OPTS:'any']
        def task = Spy(TaskRun);
        def proc = Mock(TaskProcessor)

        when:
        def env = task.getEnvironment()
        then:
        1 * task.getProcessor() >> proc
        1 * proc.getProcessEnvironment() >> [FOO: 'hola', BAR: 'world']
        1 * task.getInputEnvironment() >> [BAR: 'mundo', OMEGA: 'ooo', _OPTS: 'any']
        env ==  EXPECT   // note: `BAR` in the process config should be overridden by `BAR` in the task input
        str(env) == str(EXPECT)
    }

    private String str(Map env) {
        def result = ''
        env.each { k, v ->
            result += "$k=$v\n"
        }
        return result
    }

    def 'should get task env as string' () {

        given:
        def task = Spy(TaskRun);

        when:
        def env = task.getEnvironmentStr()
        then:
        1 * task.getEnvironment() >> [FOO: 1, BAR: 2, BAZ: 'hello world']
        env ==  '''
                FOO=1
                BAR=2
                BAZ=hello world
                '''
                .stripIndent().leftTrim()


        when:
        env = task.getEnvironmentStr()
        then:
        1 * task.getEnvironment() >> [:]
        env == null

    }


}
