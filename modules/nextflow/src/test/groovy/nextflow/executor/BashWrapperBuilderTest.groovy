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

package nextflow.executor

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.container.DockerBuilder
import nextflow.container.SingularityBuilder
import nextflow.processor.TaskBean
import nextflow.util.MustacheTemplateEngine
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderTest extends Specification {

    def setupSpec() {
        new Session()
    }

    private String load(String name, Map<String,String> binding=[:]) {
        def template = new File("src/test/groovy/nextflow/executor/$name").text
        return binding ? new MustacheTemplateEngine().render(template, binding) : template
    }

    private BashWrapperBuilder newBashWrapperBuilder(Map bean=[:]) {
        if( !bean.containsKey('workDir') )
            bean.workDir = Paths.get('/work/dir')
        if( !bean.script )
            bean.script = 'echo Hello world!'
        new BashWrapperBuilder(bean as TaskBean)
    }

    def 'test map constructor'() {

        given:
        def bean = new TaskBean(
                input: 'alpha',
                scratch: '$var_x',
                workDir: Paths.get('a'),
                targetDir: Paths.get('b'),
                containerImage: 'docker_x',
                environment: [a:1, b:2],
                script: 'echo ciao',
                shell: ['bash','-e']
        )

        when:
        def wrapper = new BashWrapperBuilder(bean)

        then:
        wrapper.scratch == '$var_x'
        wrapper.input == 'alpha'
        wrapper.workDir == Paths.get('a')
        wrapper.targetDir == Paths.get('b')
        wrapper.containerImage == 'docker_x'
        wrapper.environment ==  [a:1, b:2]
        wrapper.script ==  'echo ciao'
        wrapper.shell == ['bash','-e']
    }


    def 'should create bash launcher files' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        newBashWrapperBuilder(
                workDir: folder,
                script: 'echo Hello world!' ).build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        !Files.exists(folder.resolve('.command.in'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                ''' .stripIndent().leftTrim()

        folder.resolve('.command.run').text.contains('nxf_main')
        !folder.resolve('.command.run').text.contains('nxf_trace')

        cleanup:
        folder?.deleteDir()
    }

    def 'should create bash launcher files with trace' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        newBashWrapperBuilder(
                workDir: folder,
                statsEnabled: true ) .build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        !Files.exists(folder.resolve('.command.in'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''.stripIndent().leftTrim()

        folder.resolve('.command.run').text.contains('nxf_main')
        folder.resolve('.command.run').text.contains('nxf_trace')

        cleanup:
        folder?.deleteDir()
    }

    def 'should create launcher with input' () {
        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        newBashWrapperBuilder(
                workDir: folder,
                input: 'foo bar' ) .build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        Files.exists(folder.resolve('.command.in'))

        folder.resolve('.command.in').text == 'foo bar'
        folder.resolve('.command.sh').text.contains('echo Hello world!')
        folder.resolve('.command.run').text.contains('nxf_main')

        cleanup:
        folder?.deleteDir()
    }


    def 'should create module command' () {
        given:
        def wrapper = new BashWrapperBuilder(Mock(TaskBean))

        expect:
        wrapper.moduleLoad('foo')  == 'nxf_module_load foo'
        wrapper.moduleLoad('foo/1.2')  == 'nxf_module_load foo 1.2'
        wrapper.moduleLoad('foo/bar/1.2')  == 'nxf_module_load foo/bar 1.2'
        wrapper.moduleLoad('foo/bar/')  == 'nxf_module_load foo/bar '
    }


    def 'should create container env' () {
        given:
        def bash = Spy(BashWrapperBuilder)

        when:
        def builder = bash.createContainerBuilder(null)
        then:
        bash.getEnvironment() >> [:]
        bash.getBinDir() >> Paths.get('/my/bin')
        bash.getWorkDir() >> Paths.get('/my/work/dir')
        bash.getStatsEnabled() >> false
        bash.getStageInMode() >> 'symlink'

        bash.getInputFiles() >> [:]
        bash.getContainerConfig() >> [engine: 'singularity', envWhitelist: 'FOO,BAR']
        bash.getContainerImage() >> 'foo/bar'
        bash.getContainerMount() >> null
        bash.getContainerMemory() >> null
        bash.getContainerCpuset() >> null
        bash.getContainerOptions() >> null

        builder instanceof SingularityBuilder
        builder.env == ['FOO','BAR']
        builder.workDir == Paths.get('/my/work/dir')
    }

    def 'should add resolved inputs'() {
        given:
        def bash = Spy(BashWrapperBuilder)
        bash.bean = Mock(TaskBean)
        bash.getContainerConfig() >> [engine: 'docker']

        def BUILDER = Mock(DockerBuilder)
        def INPUTS = ['/some/path': Paths.get('/store/path.txt')]

        // check input files are mounted in the container
        when:
        bash.createContainerBuilder(null)
        then:
        bash.createContainerBuilder0('docker') >> BUILDER
        bash.getInputFiles() >> INPUTS
        bash.getStageInMode() >> null
        1 * BUILDER.addMountForInputs(INPUTS) >> null

        // -- do not mount inputs when stage-in-mode == 'copy'
        when:
        bash.createContainerBuilder(null)
        then:
        bash.createContainerBuilder0('docker') >> BUILDER
        bash.getStageInMode() >> 'copy'
        0 * BUILDER.addMountForInputs(_) >> null
    }

    def 'should render launcher script' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def wrapper = newBashWrapperBuilder(
                name: 'Hello 1',
                workDir: folder,
                headerScript: '#BSUB -x 1\n#BSUB -y 2',
                beforeScript: 'echo Before',
                condaEnv: Paths.get('/conda/env/path') ) .buildNew0()

        then:
        wrapper == load('test-bash-wrapper.txt', [folder: folder.toString()])

        cleanup:
        folder?.deleteDir()
    }

    def 'should render launcher script with trace' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def wrapper = newBashWrapperBuilder(
                name: 'Hello 2',
                workDir: folder,
                statsEnabled: true ) .buildNew0()

        then:
        wrapper == load('test-bash-wrapper-with-trace.txt', [folder: folder.toString()])

        cleanup:
        folder?.deleteDir()
    }

    @Unroll
    def 'test change to scratchDir' () {

        setup:
        def builder = newBashWrapperBuilder()

        when:
        builder.scratch = SCRATCH
        then:
        builder.makeBinding().scratch_cmd == EXPECTED


        where:
        SCRATCH     | EXPECTED
        null        | "NXF_SCRATCH=''"
        true        | 'NXF_SCRATCH="$(set +u; nxf_mktemp $TMPDIR)"'
        '$SOME_DIR' | 'NXF_SCRATCH="$(set +u; nxf_mktemp $SOME_DIR)"'
        '/my/temp'  | 'NXF_SCRATCH="$(set +u; nxf_mktemp /my/temp)"'
        'ram-disk'  | 'NXF_SCRATCH="$(nxf_mktemp /dev/shm)"'

    }

    def 'should return task name' () {
        expect:
        newBashWrapperBuilder(name: 'foo').makeBinding().task_name == 'foo'
        newBashWrapperBuilder(name: 'bar').makeBinding().task_name == 'bar'
    }

    def 'should return header directives' () {
        when:
        def bash = newBashWrapperBuilder()
        then:
        bash.makeBinding().containsKey('header_script')
        bash.makeBinding().header_script == null

        when:
        bash = newBashWrapperBuilder(headerScript: '#BSUB -x 1\n#BSUB -y 2')
        then:
        bash.makeBinding().header_script == '#BSUB -x 1\n#BSUB -y 2'

    }

    def 'should copy control files' () {

        when:
        def binding = newBashWrapperBuilder(scratch: false).makeBinding()
        then:
        binding.containsKey('unstage_controls')
        binding.unstage_controls == null

        when:
        binding = newBashWrapperBuilder(scratch: true).makeBinding()
        then:
        binding.unstage_controls == '''\
                cp .command.out /work/dir/.command.out || true
                cp .command.err /work/dir/.command.err || true
                '''.stripIndent()

        when:
        binding = newBashWrapperBuilder(scratch: true, statsEnabled: true).makeBinding()
        then:
        binding.unstage_controls == '''\
                cp .command.out /work/dir/.command.out || true
                cp .command.err /work/dir/.command.err || true
                cp .command.trace /work/dir/.command.trace || true
                '''.stripIndent()

    }

    def 'should stage inputs' () {

        given:
        def folder = Paths.get('/work/dir')
        def inputs = ['sample_1.fq':Paths.get('/some/data/sample_1.fq'), 'sample_2.fq':Paths.get('/some/data/sample_2.fq'), ]

        when:
        def binding = newBashWrapperBuilder([
                workDir: folder,
                targetDir: folder,
                inputFiles: inputs ]).makeBinding()

        then:
        binding.stage_inputs == '''\
                # stage input files
                rm -f sample_1.fq
                rm -f sample_2.fq
                ln -s /some/data/sample_1.fq sample_1.fq
                ln -s /some/data/sample_2.fq sample_2.fq
                '''.stripIndent().rightTrim()


    }

    def 'should unstage outputs' () {

        given:
        def folder = Paths.get('/work/dir')
        def outputs = ['test.bam','test.bai']

        when:
        def binding = newBashWrapperBuilder([
                workDir: folder,
                targetDir: folder,
                scratch: false,
                outputFiles: outputs ]).makeBinding()

        then:
        binding.containsKey('unstage_outputs')
        binding.unstage_outputs == null


        when:
        binding = newBashWrapperBuilder([
                workDir: folder,
                targetDir: folder,
                scratch: true,
                outputFiles: outputs ]).makeBinding()

        then:
        binding.unstage_outputs == '''\
                mkdir -p /work/dir
                cp -fRL test.bam /work/dir || true
                cp -fRL test.bai /work/dir || true
                '''.stripIndent().rightTrim()


        when:
        binding = newBashWrapperBuilder([
                workDir: folder,
                targetDir: Paths.get('/another/dir'),
                scratch: false,
                outputFiles: outputs ]).makeBinding()

        then:
        binding.unstage_outputs == '''\
                mkdir -p /another/dir
                mv -f test.bam /another/dir || true
                mv -f test.bai /another/dir || true
                '''.stripIndent().rightTrim()
    }

    def 'should create env' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.containsKey('task_env')
        binding.containsKey('container_env')
        binding.task_env == null
        binding.container_env == null

        when:
        binding = newBashWrapperBuilder(environment: [FOO:'aa', BAR:'bb']).makeBinding()
        then:
        binding.containsKey('task_env')
        binding.containsKey('container_env')
        binding.task_env == /\
                export FOO="aa"
                export BAR="bb"
                /.stripIndent()
        binding.container_env == null

        when:
        binding = newBashWrapperBuilder(environment: [FOO:'aa', BAR:'bb'], containerEnabled: true, containerImage: 'foo', containerConfig: [engine: 'docker']).makeBinding()
        then:
        binding.containsKey('task_env')
        binding.containsKey('container_env')
        binding.task_env == null
        binding.container_env == /\
                nxf_container_env() {
                cat << EOF
                export FOO="aa"
                export BAR="bb"
                EOF
                }
                /.stripIndent()
        binding.launch_cmd.contains('nxf_container_env')

        when:
        binding = newBashWrapperBuilder(environment: [FOO:'aa', BAR:'bb'], containerEnabled: true, containerImage: 'foo', containerNative: true).makeBinding()
        then:
        binding.containsKey('task_env')
        binding.containsKey('container_env')
        binding.task_env == 'export FOO="aa"\nexport BAR="bb"\n'
        binding.container_env == null
        !binding.launch_cmd.contains('nxf_container_env')
    }

    def 'should enable trace feature' () {
        when:
        def binding = newBashWrapperBuilder(statsEnabled: false).makeBinding()
        then:
        binding.launch_cmd == '/bin/bash -ue /work/dir/.command.sh'
        binding.unstage_controls == null
        binding.containsKey('unstage_controls')

        when:
        binding = newBashWrapperBuilder(statsEnabled: true).makeBinding()
        then:
        binding.launch_cmd == '/bin/bash /work/dir/.command.run nxf_trace'
        binding.unstage_controls == null
        binding.containsKey('unstage_controls')

        when:
        binding = newBashWrapperBuilder(statsEnabled: true, scratch: true).makeBinding()
        then:
        binding.launch_cmd == '/bin/bash /work/dir/.command.run nxf_trace'
        binding.unstage_controls == '''\
                        cp .command.out /work/dir/.command.out || true
                        cp .command.err /work/dir/.command.err || true
                        cp .command.trace /work/dir/.command.trace || true
                        '''.stripIndent()

    }

    def 'should create launcher command with input' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.launch_cmd == '/bin/bash -ue /work/dir/.command.sh'
        binding.trace_cmd == binding.launch_cmd

        when:
        binding = newBashWrapperBuilder(input: 'Ciao ciao').makeBinding()
        then:
        binding.launch_cmd == '/bin/bash -ue /work/dir/.command.sh < /work/dir/.command.in'
        binding.trace_cmd == binding.launch_cmd

        when:
        def PERL = '''
            #!/usr/bin/env perl
            print "Hello world\n";
            '''
        binding = newBashWrapperBuilder(input: 'Ciao ciao', script: PERL).makeBinding()
        then:
        binding.launch_cmd == '/usr/bin/env perl /work/dir/.command.sh < /work/dir/.command.in'
        binding.trace_cmd == binding.launch_cmd
    }

    def 'should create before script' () {
        when:
        def binding = newBashWrapperBuilder(beforeScript: 'do_this_and_that').makeBinding()
        then:
        binding.before_script == '''\
                # beforeScript directive
                do_this_and_that
                '''.stripIndent()

        when:
        binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.before_script == null
        binding.containsKey('before_script')
    }

    def 'should create module load snippet' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.module_load == null

        when:
        binding = newBashWrapperBuilder(moduleNames: ['ciao/1', 'mondo/2', 'bioinfo-tools']).makeBinding()
        then:
        binding.module_load == '''\
                nxf_module_load ciao 1
                nxf_module_load mondo 2
                nxf_module_load bioinfo-tools
                '''.stripIndent()

        binding.helpers_script == '''\
                nxf_module_load(){
                  local mod=$1
                  local ver=${2:-}
                  local new_module="$mod"; [[ $ver ]] && new_module+="/$ver"
                
                  if [[ ! $(module list 2>&1 | grep -o "$new_module") ]]; then
                    old_module=$(module list 2>&1 | grep -Eow "$mod\\/[^\\( \\n]+" || true)
                    if [[ $ver && $old_module ]]; then
                      module switch $old_module $new_module
                    else
                      module load $new_module
                    fi
                  fi
                }
                
                '''.stripIndent()

    }

    def 'should create conda activate snippet' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.conda_activate == null
        binding.containsKey('conda_activate')

        when:
        def CONDA = Paths.get('/some/conda/env/foo')
        binding = newBashWrapperBuilder(condaEnv: CONDA).makeBinding()
        then:
        binding.conda_activate == '''\
                # conda environment
                source $(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }')/bin/activate /some/conda/env/foo
                '''.stripIndent()

    }

    def 'should cleanup scratch dir' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.cleanup_cmd == ''
        binding.containsKey('cleanup_cmd')


        when:
        binding = newBashWrapperBuilder(scratch: true).makeBinding()
        then:
        binding.cleanup_cmd == 'rm -rf $NXF_SCRATCH || true\n'
    }

    def 'should create wrapper with docker' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerEnabled: true,
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash --name $NXF_BOXID busybox -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker kill $NXF_BOXID'
    }

    def 'should create wrapper with docker with sudo' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerConfig: [engine: 'docker', sudo: true, enabled: true],
                containerEnabled: true  ).makeBinding()

        then:
        binding.launch_cmd == 'sudo docker run -i -v /work/dir:/work/dir -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash --name $NXF_BOXID busybox -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == 'sudo docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'sudo docker kill $NXF_BOXID'
    }

    def 'should create wrapper with docker no kill' () {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'ubuntu',
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true, remove:false, kill: false] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash --name $NXF_BOXID ubuntu -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == null
        binding.containsKey('kill_cmd')
    }

    def 'should create wrapper with docker with custom signal' () {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'ubuntu',
                containerConfig: [engine: 'docker', enabled: true, remove:false, kill: 'SIGXXX'] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -v /work/dir:/work/dir -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash --name $NXF_BOXID ubuntu -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == 'docker kill -s SIGXXX $NXF_BOXID'
        binding.containsKey('kill_cmd')
    }

    def 'should create wrapper with docker with mount' () {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'busybox',
                containerMount: '/folder with blanks' as Path,
                containerConfig: [engine: 'docker', enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v /work/dir:/work/dir -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID busybox -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker kill $NXF_BOXID'
    }

    def 'should create wrapper with docker with scratch' () {
        when:
        def binding = newBashWrapperBuilder(
                scratch: true,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'busybox',
                containerConfig: [engine: 'docker', sudo: true, enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'sudo docker run -i -v /work/dir:/work/dir -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash --name $NXF_BOXID busybox -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.kill_cmd == 'sudo docker kill $NXF_BOXID'
        binding.cleanup_cmd == '''\
                (sudo -n true && sudo rm -rf "$NXF_SCRATCH" || rm -rf "$NXF_SCRATCH")&>/dev/null || true
                sudo docker rm $NXF_BOXID &>/dev/null || true
                '''.stripIndent()
    }

    def 'should create wrapper with docker container custom options' () {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'busybox',
                containerOptions: '-v /foo:/bar',
                containerConfig: [engine: 'docker', enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -v /work/dir:/work/dir -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash -v /foo:/bar --name $NXF_BOXID busybox -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.kill_cmd == 'docker kill $NXF_BOXID'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
    }

    def 'should create wrapper with shifter'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'docker:ubuntu:latest',
                environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'shifter'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == '''\
        shifterimg pull docker:ubuntu:latest
        shifterimg lookup docker:ubuntu:latest
        while ! shifterimg lookup docker:ubuntu:latest; do
            sleep 5
            STATUS=$(shifterimg -v pull docker:ubuntu:latest | tail -n2 | head -n1 | awk \'{print $6}\')
            [[ $STATUS == "FAILURE" || -z $STATUS ]] && echo "Shifter failed to pull image \'docker:ubuntu:latest\'" >&2  && exit 1
        done
        shifter --image docker:ubuntu:latest /bin/bash -c "eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"
        '''.stripIndent().rightTrim()
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && kill $pid 2>/dev/null'

    }

    def 'should create wrapper with singularity'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'docker:ubuntu:latest',
                environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'singularity'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec docker:ubuntu:latest /bin/bash -c "cd $PWD; eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && kill $pid 2>/dev/null'

    }

    def 'should create task and container env' () {
        given:
        def ENV = [FOO: 'hello', BAR: 'hello world', PATH: '/some/path:$PATH']

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.task_env == null
        binding.container_env == null

        when:
        binding = newBashWrapperBuilder(environment: ENV).makeBinding()
        then:
        binding.container_env == null
        binding.task_env ==  '''
                export FOO="hello"
                export BAR="hello world"
                export PATH="/some/path:$PATH"
                '''
                .stripIndent().leftTrim()

        when:
        binding = newBashWrapperBuilder(environment: ENV,
                containerEnabled: true,
                containerImage: 'busybox',
                containerConfig: [enabled: true, engine: 'docker']).makeBinding()
        then:
        binding.task_env == null
        binding.container_env ==  '''
                nxf_container_env() {
                cat << EOF
                export FOO="hello"
                export BAR="hello world"
                export PATH="/some/path:\\$PATH"
                EOF
                }
                '''
                .stripIndent().leftTrim()
    }

    def 'should include fix ownership command' () {

        given:
        def bean = Mock(TaskBean)
        def copy = Mock(ScriptFileCopyStrategy)
        bean.workDir >> Paths.get('/work/dir')
        def builder = Spy(BashWrapperBuilder)
        builder.bean = bean
        builder.copyStrategy = copy

        when:
        def binding = builder.makeBinding()
        then:
        builder.fixOwnership() >> true
        binding.fix_ownership == '[ ${NXF_OWNER:=\'\'} ] && chown -fR --from root $NXF_OWNER /work/dir/{*,.*} || true'


        when:
        binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.fix_ownership == null
        binding.containsKey('fix_ownership')
    }

    def 'should create script command' () {
        when:
        def binding = newBashWrapperBuilder(afterScript: "cleanup_that" ).makeBinding()
        then:
        binding.after_script == "# 'afterScript' directive\ncleanup_that"

        when:
        binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.after_script == null
        binding.containsKey('after_script')
    }


    def 'should get output env capture snippet' () {
        given:
        def builder = new BashWrapperBuilder()

        when:
        def str = builder.getOutputEnvCaptureSnippet(['FOO','BAR'])
        then:
        str == '''
            # capture process environment
            set +u
            echo FOO=$FOO > .command.env
            echo BAR=$BAR >> .command.env
            '''
            .stripIndent()

    }

    def 'should validate bash interpreter' () {
        given:
        def builder = new BashWrapperBuilder()
        expect:
        builder.isBash('/bin/bash')
        builder.isBash('/usr/bin/bash')
        builder.isBash('/bin/env bash')
        builder.isBash('/bin/bash -eu')
        !builder.isBash('/bin/env perl')

    }

    def 'should get stage and unstage commands' () {

        when:
        def builder = newBashWrapperBuilder()
        then:
        builder.getStageCommand() == 'nxf_stage'
        builder.getUnstageCommand() == 'nxf_unstage'

        when:
        def binding = builder.makeBinding()
        then:
        binding.stage_cmd == 'nxf_stage'
        binding.unstage_cmd == 'nxf_unstage'
        
    }
}
