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

package nextflow.executor

import java.nio.file.FileSystemException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Session
import nextflow.SysEnv
import nextflow.container.ContainerConfig
import nextflow.container.DockerBuilder
import nextflow.container.SingularityBuilder
import nextflow.processor.TaskBean
import nextflow.util.MustacheTemplateEngine
import org.yaml.snakeyaml.Yaml
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
        def template = new File("src/test/resources/nextflow/executor/$name").text
        return binding ? new MustacheTemplateEngine().render(template, binding) : template
    }

    private BashWrapperBuilder newBashWrapperBuilder(Map bean=[:]) {
        if( !bean.containsKey('workDir') )
            bean.workDir = Paths.get('/work/dir')
        if( !bean.script )
            bean.script = 'echo Hello world!'
        if( !bean.containsKey('inputFiles') )
            bean.inputFiles = [:]
        if( !bean.containsKey('outputFiles') )
            bean.outputFiles = []
        new BashWrapperBuilder(bean as TaskBean) {
            @Override
            protected String getSecretsEnv() {
                return null
            }
        }
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
        and:
        def builder = newBashWrapperBuilder( workDir: folder, input: 'foo bar' )

        when:
        builder.build()
        then:
        builder.targetInputFile() == folder.resolve('.command.in')
        builder.targetScriptFile() == folder.resolve('.command.sh')
        builder.targetWrapperFile() == folder.resolve('.command.run')
        builder.targetStageFile() == folder.resolve('.command.stage')
        and:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        Files.exists(folder.resolve('.command.in'))
        and:
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
        and:
        bash.getEnvironment() >> [:]
        bash.getBinDirs() >> [Paths.get('/my/bin') ]
        bash.getWorkDir() >> Paths.get('/my/work/dir')
        bash.isStatsEnabled() >> false
        bash.getStageInMode() >> 'symlink'
        bash.getInputFiles() >> [:]
        bash.getContainerConfig() >> [engine: 'singularity', envWhitelist: 'FOO,BAR']
        bash.getContainerImage() >> 'foo/bar'
        bash.getContainerMount() >> null
        bash.getContainerMemory() >> null
        bash.getContainerCpus() >> null
        bash.getContainerCpuset() >> null
        bash.getContainerOptions() >> null
        bash.isSecretNative() >> false
        bash.getSecretNames() >> []

        when:
        def builder = bash.createContainerBuilder(null)
        then:
        builder instanceof SingularityBuilder
        builder.env == ['NXF_TASK_WORKDIR', 'FOO','BAR']
        builder.workDir == Paths.get('/my/work/dir')
        builder.mounts == [ Paths.get('/my/bin') ]
        
    }

    def 'should add resolved inputs'() {
        given:
        def bash = Spy(new BashWrapperBuilder(bean: Mock(TaskBean)))
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
                condaEnv: Paths.get('/conda/env/path'),
                spackEnv: Paths.get('/spack/env/path') ) .buildNew0()

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

    def 'should create task metadata string' () {
        given:
        def builder = newBashWrapperBuilder(
            name: 'foo',
            arrayIndexName: 'SLURM_ARRAY_TASK_ID',
            arrayIndexStart: 0,
            arrayWorkDirs: [ Path.of('/work/01'), Path.of('/work/02'), Path.of('/work/03') ],
            containerConfig: [enabled: true],
            containerImage: 'quay.io/nextflow:bash',
            outputFiles: ['foo.txt', '*.bar', '**/baz']
        )

        when:
        def meta = builder.getTaskMetadata()
        then:
        meta == '''\
            ### ---
            ### name: 'foo'
            ### array:
            ###   index-name: SLURM_ARRAY_TASK_ID
            ###   index-start: 0
            ###   work-dirs:
            ###   - /work/01
            ###   - /work/02
            ###   - /work/03
            ### container: 'quay.io/nextflow:bash'
            ### outputs:
            ### - 'foo.txt'
            ### - '*.bar'
            ### - '**/baz'
            ### ...
            '''.stripIndent()

        when:
        def yaml = meta.readLines().collect(it-> it.substring(4)).join('\n')
        def obj = new Yaml().load(yaml) as Map
        then:
        obj.name == 'foo'
        obj.array == [
            'index-name':'SLURM_ARRAY_TASK_ID',
            'index-start':0,
            'work-dirs':['/work/01', '/work/02', '/work/03']
        ]
        obj.container == 'quay.io/nextflow:bash'
        obj.outputs == ['foo.txt', '*.bar', '**/baz']
    }

    def 'should add task metadata' () {
        when:
        def bash = newBashWrapperBuilder([name:'task1'])
        then:
        bash.makeBinding().containsKey('task_metadata')
        bash.makeBinding().task_metadata == '''\
            ### ---
            ### name: 'task1'
            ### ...
            '''.stripIndent()

        when:
        bash = newBashWrapperBuilder(
            name: 'task2',
            arrayIndexName: 'SLURM_ARRAY_TASK_ID',
            arrayIndexStart: 0,
            arrayWorkDirs: [ Path.of('/work/01'), Path.of('/work/02'), Path.of('/work/03') ],
            containerConfig: [enabled: true],
            containerImage: 'quay.io/nextflow:bash',
            outputFiles: ['foo.txt', '*.bar', '**/baz']
        )
        then:
        bash.makeBinding().task_metadata == '''\
            ### ---
            ### name: 'task2'
            ### array:
            ###   index-name: SLURM_ARRAY_TASK_ID
            ###   index-start: 0
            ###   work-dirs:
            ###   - /work/01
            ###   - /work/02
            ###   - /work/03
            ### container: 'quay.io/nextflow:bash'
            ### outputs:
            ### - 'foo.txt'
            ### - '*.bar'
            ### - '**/baz'
            ### ...
            '''.stripIndent()
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
        def inputs = [
            'sample_1.fq': Paths.get('/some/data/sample_1.fq'),
            'sample_2.fq': Paths.get('/some/data/sample_2.fq'),
        ]
        def stageScript = '''\
                # stage input files
                rm -f sample_1.fq
                rm -f sample_2.fq
                ln -s /some/data/sample_1.fq sample_1.fq
                ln -s /some/data/sample_2.fq sample_2.fq
                '''.stripIndent().rightTrim()

        when:
        def binding = newBashWrapperBuilder([
                workDir: folder,
                targetDir: folder,
                inputFiles: inputs ]).makeBinding()

        then:
        binding.stage_inputs == stageScript
    }

    def 'should stage inputs to external file' () {
        given:
        SysEnv.push([NXF_WRAPPER_STAGE_FILE_THRESHOLD: '100'])
        and:
        def folder = Files.createTempDirectory('test')
        and:
        def inputs = [
                'sample_1.fq': Paths.get('/some/data/sample_1.fq'),
                'sample_2.fq': Paths.get('/some/data/sample_2.fq'),
        ]
        def stageScript = '''\
                rm -f sample_1.fq
                rm -f sample_2.fq
                ln -s /some/data/sample_1.fq sample_1.fq
                ln -s /some/data/sample_2.fq sample_2.fq
                '''.stripIndent().rightTrim()
        and:
        def builder = newBashWrapperBuilder([
                workDir: folder,
                targetDir: folder,
                inputFiles: inputs ])

        when:
        def binding = builder.makeBinding()
        then:
        binding.stage_inputs == "# stage input files\nbash ${folder}/.command.stage"

        when:
        builder.build()
        then:
        folder.resolve('.command.stage').text == stageScript

        cleanup:
        SysEnv.pop()
        folder?.deleteDir()
    }

    def 'should include sync command' () {
        given:
        SysEnv.push([NXF_ENABLE_FS_SYNC: 'true'])

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.sync_cmd == 'sync || true'

        cleanup:
        SysEnv.pop()
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
                IFS=$'\\n'
                for name in $(eval "ls -1d test.bam test.bai" | sort | uniq); do
                    nxf_fs_copy "$name" /work/dir
                done
                unset IFS
                '''.stripIndent().rightTrim()


        when:
        binding = newBashWrapperBuilder([
                workDir: folder,
                targetDir: Paths.get('/another/dir'),
                scratch: false,
                outputFiles: outputs ]).makeBinding()

        then:
        binding.unstage_outputs == '''\
                IFS=$'\\n'
                for name in $(eval "ls -1d test.bam test.bai" | sort | uniq); do
                    nxf_fs_move "$name" /another/dir
                done
                unset IFS
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
        binding.sync_cmd == null

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
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create secret env' () {

        when:
        def binding0 = newBashWrapperBuilder().makeBinding()
        then:
        binding0.containsKey('secrets_env')
        binding0.secrets_env == null

        when:
        def builder1 = Spy(newBashWrapperBuilder()) {
            getSecretsEnv() >> 'source /some/file.txt'
            isSecretNative() >> false
        }
        and:
        def binding1 = builder1.makeBinding()
        then:
        binding1.secrets_env == 'source /some/file.txt'

        when:
        def builder2 = Spy(newBashWrapperBuilder()) {
            getSecretsEnv() >> 'source /some/file.txt'
            isSecretNative() >> true
        }
        and:
        def binding2 = builder2.makeBinding()
        then:
        binding2.secrets_env == null
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

    def 'should enable trace feature with custom shell' () {
        when:
        def binding = newBashWrapperBuilder(statsEnabled: false, shell: ['/bin/bash', '-ue']).makeBinding()
        then:
        binding.launch_cmd == '/bin/bash -ue /work/dir/.command.sh'
        binding.unstage_controls == null
        binding.containsKey('unstage_controls')

        when:
        binding = newBashWrapperBuilder(statsEnabled: true, shell: ['/bin/bash', '-eu']).makeBinding()
        then:
        binding.launch_cmd == '/bin/bash -eu /work/dir/.command.run nxf_trace'
        binding.unstage_controls == null
        binding.containsKey('unstage_controls')

        when:
        binding = newBashWrapperBuilder(statsEnabled: true, scratch: true, shell: ['/bin/bash', '-eu']).makeBinding()
        then:
        binding.launch_cmd == '/bin/bash -eu /work/dir/.command.run nxf_trace'
        binding.unstage_controls == '''\
                        cp .command.out /work/dir/.command.out || true
                        cp .command.err /work/dir/.command.err || true
                        cp .command.trace /work/dir/.command.trace || true
                        '''.stripIndent()

        when:
        binding = newBashWrapperBuilder(statsEnabled: true, shell: ['/usr/local/bin/bash', '-ue']).makeBinding()
        then:
        binding.launch_cmd == '/usr/local/bin/bash -ue /work/dir/.command.run nxf_trace'

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

    def 'should create micromamba activate snippet' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.conda_activate == null
        binding.containsKey('conda_activate')

        when:
        def CONDA = Paths.get('/some/conda/env/foo')
        binding = newBashWrapperBuilder([condaEnv: CONDA, 'useMicromamba': true]).makeBinding()
        then:
        binding.conda_activate == '''\
                # conda environment
                eval "$(micromamba shell hook --shell bash)" && micromamba activate /some/conda/env/foo
                '''.stripIndent()
    }

    def 'should create spack activate snippet' () {

        when:
        def binding = newBashWrapperBuilder().makeBinding()
        then:
        binding.spack_activate == null
        binding.containsKey('spack_activate')

        when:
        def SPACK = Paths.get('/some/spack/env/foo')
        binding = newBashWrapperBuilder(spackEnv: SPACK).makeBinding()
        then:
        binding.spack_activate == '''\
                # spack environment
                spack env activate -d /some/spack/env/foo
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
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker stop $NXF_BOXID'
    }

    def 'should create wrapper with docker and environment' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerEnabled: true,
                environment: [FOO: 'something'],
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -c "eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker stop $NXF_BOXID'
        and:
        binding.container_env == '''\
                    nxf_container_env() {
                    cat << EOF
                    export FOO="something"
                    EOF
                    }
                    '''.stripIndent()
    }

    def 'should create wrapper with docker and legacy entrypoint false' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerEnabled: true,
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true, entrypointOverride: false] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker stop $NXF_BOXID'
    }

    def 'should create wrapper with docker and legacy entrypoint false and env' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerEnabled: true,
                environment: [FOO:'hello'],
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true, entrypointOverride: false] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -c "eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker stop $NXF_BOXID'
        and:
        binding.container_env == '''\
                    nxf_container_env() {
                    cat << EOF
                    export FOO="hello"
                    EOF
                    }
                    '''.stripIndent()
    }

    def 'should create wrapper with docker with sudo' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerConfig: [engine: 'docker', sudo: true, enabled: true],
                containerEnabled: true  ).makeBinding()

        then:
        binding.launch_cmd == 'sudo docker run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.cleanup_cmd == 'sudo docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'sudo docker stop $NXF_BOXID'
    }

    def 'should create wrapper with docker no kill' () {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'ubuntu',
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true, remove:false, kill: false] ).makeBinding()

        then:
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v $(nxf_mktemp):/tmp -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID ubuntu /bin/bash -ue /work/dir/.command.sh'
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
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID ubuntu /bin/bash -ue /work/dir/.command.sh'
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
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v /work/dir:/work/dir -w "\$NXF_TASK_WORKDIR" --name \$NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'docker stop $NXF_BOXID'
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
        binding.launch_cmd == 'sudo docker run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.kill_cmd == 'sudo docker stop $NXF_BOXID'
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
        binding.launch_cmd == 'docker run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" -v /foo:/bar --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.kill_cmd == 'docker stop $NXF_BOXID'
        binding.cleanup_cmd == 'docker rm $NXF_BOXID &>/dev/null || true\n'
    }

    def 'should create wrapper with sarus'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'busybox',
                containerConfig: [enabled: true, engine: 'sarus'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == '''\
        sarus pull busybox 1>&2
        sarus run -e "NXF_TASK_WORKDIR" --mount=type=bind,source=/work/dir,destination=/work/dir -w "$NXF_TASK_WORKDIR" busybox /bin/bash -ue /work/dir/.command.sh
        '''.stripIndent().rightTrim()
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create wrapper with sarus and environment'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'busybox',
                environment: [FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'sarus'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == '''\
        sarus pull busybox 1>&2
        sarus run -e "NXF_TASK_WORKDIR" --mount=type=bind,source=/work/dir,destination=/work/dir -w "$NXF_TASK_WORKDIR" busybox /bin/bash -c "eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"
        '''.stripIndent().rightTrim()
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
        and:
        binding.container_env == '''\
                    nxf_container_env() {
                    cat << EOF
                    export FOO="xxx"
                    EOF
                    }
                    '''.stripIndent()
    }

    def 'should create wrapper with sarus with mount'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'busybox',
                containerMount: '/folder with blanks' as Path,
                containerConfig: [enabled: true, engine: 'sarus'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == '''\
        sarus pull busybox 1>&2
        sarus run -e "NXF_TASK_WORKDIR" --mount=type=bind,source=/folder\\ with\\ blanks,destination=/folder\\ with\\ blanks --mount=type=bind,source=/work/dir,destination=/work/dir -w "$NXF_TASK_WORKDIR" busybox /bin/bash -ue /work/dir/.command.sh
        '''.stripIndent().rightTrim()
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create wrapper with sarus container custom options'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'busybox',
                containerOptions: '--mount=type=bind,source=/foo,destination=/bar',
                containerConfig: [enabled: true, engine: 'sarus'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == '''\
        sarus pull busybox 1>&2
        sarus run -e "NXF_TASK_WORKDIR" --mount=type=bind,source=/work/dir,destination=/work/dir -w "$NXF_TASK_WORKDIR" --mount=type=bind,source=/foo,destination=/bar busybox /bin/bash -ue /work/dir/.command.sh
        '''.stripIndent().rightTrim()
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create wrapper with shifter'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'docker://ubuntu:latest',
                environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'shifter'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == '''\
        shifterimg pull docker://ubuntu:latest
        shifterimg lookup docker://ubuntu:latest
        while ! shifterimg lookup docker://ubuntu:latest; do
            sleep 5
            STATUS=$(shifterimg -v pull docker://ubuntu:latest | tail -n2 | head -n1 | awk \'{print $6}\')
            [[ $STATUS == "FAILURE" || -z $STATUS ]] && echo "Shifter failed to pull image \'docker://ubuntu:latest\'" >&2  && exit 1
        done
        ${NXF_TASK_WORKDIR:+"NXF_TASK_WORKDIR=$NXF_TASK_WORKDIR"} shifter --image docker://ubuntu:latest /bin/bash -c "eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"
        '''.stripIndent().rightTrim()
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'

    }

    def 'should create wrapper with singularity'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'docker://ubuntu:latest',
                environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'singularity'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${NXF_TASK_WORKDIR:+SINGULARITYENV_NXF_TASK_WORKDIR="$NXF_TASK_WORKDIR"} singularity exec --no-home --pid -B /work/dir docker://ubuntu:latest /bin/bash -c "cd $NXF_TASK_WORKDIR; eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create wrapper with singularity and no env'() {
        when:
        def binding = newBashWrapperBuilder(
            containerEnabled: true,
            containerImage: 'docker://ubuntu:latest',
            environment: [:],
            containerConfig: [enabled: true, engine: 'singularity'] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${NXF_TASK_WORKDIR:+SINGULARITYENV_NXF_TASK_WORKDIR="$NXF_TASK_WORKDIR"} singularity exec --no-home --pid -B /work/dir docker://ubuntu:latest /bin/bash -c "cd $NXF_TASK_WORKDIR; /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create wrapper with singularity legacy entry'() {
        when:
        def binding = newBashWrapperBuilder(
                containerEnabled: true,
                containerImage: 'docker://ubuntu:latest',
                environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'singularity', entrypointOverride: true] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${NXF_TASK_WORKDIR:+SINGULARITYENV_NXF_TASK_WORKDIR="$NXF_TASK_WORKDIR"} singularity exec --no-home --pid -B /work/dir docker://ubuntu:latest /bin/bash -c "cd $NXF_TASK_WORKDIR; eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
    }

    def 'should create wrapper with singularity oci mode'() {
        when:
        def binding = newBashWrapperBuilder(
            containerEnabled: true,
            containerImage: 'docker://ubuntu:latest',
            environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
            containerConfig: [enabled: true, engine: 'singularity', oci: true] as ContainerConfig ).makeBinding()

        then:
        binding.launch_cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${XDG_RUNTIME_DIR:+XDG_RUNTIME_DIR="$XDG_RUNTIME_DIR"} ${DBUS_SESSION_BUS_ADDRESS:+DBUS_SESSION_BUS_ADDRESS="$DBUS_SESSION_BUS_ADDRESS"} ${NXF_TASK_WORKDIR:+SINGULARITYENV_NXF_TASK_WORKDIR="$NXF_TASK_WORKDIR"} singularity exec --no-home --oci -B /work/dir docker://ubuntu:latest /bin/bash -c "cd $NXF_TASK_WORKDIR; eval $(nxf_container_env); /bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == ""
        binding.kill_cmd == '[[ "$pid" ]] && nxf_kill $pid'
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
        def bean = Mock(TaskBean) {
            inputFiles >> [:]
            shell >> BashWrapperBuilder.BASH
            outputFiles >> []
        }
        def copy = Mock(ScriptFileCopyStrategy)
        bean.workDir >> Paths.get('/work/dir')
        and:
        def builder = Spy(new BashWrapperBuilder(bean:bean))
        builder.copyStrategy = copy

        when:
        def binding = builder.makeBinding()
        then:
        builder.fixOwnership() >> true
        binding.fix_ownership == '[ ${NXF_OWNER:=\'\'} ] && (shopt -s extglob; GLOBIGNORE=\'..\'; chown -fR --from root $NXF_OWNER /work/dir/{*,.*}) || true'


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
        def str = builder.getOutputEnvCaptureSnippet(['FOO','BAR'], Map.of())
        then:
        str == '''
            # capture process environment
            set +u
            set +e
            cd "$NXF_TASK_WORKDIR"
            
            nxf_eval_cmd() {
                {
                    IFS=$'\\n' read -r -d '' "${1}";
                    IFS=$'\\n' read -r -d '' "${2}";
                    (IFS=$'\\n' read -r -d '' _ERRNO_; return ${_ERRNO_});
                } < <((printf '\\0%s\\0%d\\0' "$(((({ shift 2; "${@}"; echo "${?}" 1>&3-; } | tr -d '\\0' 1>&4-) 4>&2- 2>&1- | tr -d '\\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
            }
            
            echo '' > .command.env
            #
            echo FOO="${FOO[@]}" >> .command.env
            echo /FOO/ >> .command.env
            #
            echo BAR="${BAR[@]}" >> .command.env
            echo /BAR/ >> .command.env
            '''
            .stripIndent()
    }

    def 'should return env & cmd capture snippet' () {
        given:
        def builder = new BashWrapperBuilder()

        when:
        def str = builder.getOutputEnvCaptureSnippet(['FOO'], [THIS: 'this --cmd', THAT: 'other "quoted" --cmd'])
        then:
        str == '''
            # capture process environment
            set +u
            set +e
            cd "$NXF_TASK_WORKDIR"
            
            nxf_eval_cmd() {
                {
                    IFS=$'\\n' read -r -d '' "${1}";
                    IFS=$'\\n' read -r -d '' "${2}";
                    (IFS=$'\\n' read -r -d '' _ERRNO_; return ${_ERRNO_});
                } < <((printf '\\0%s\\0%d\\0' "$(((({ shift 2; "${@}"; echo "${?}" 1>&3-; } | tr -d '\\0' 1>&4-) 4>&2- 2>&1- | tr -d '\\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
            }
            
            echo '' > .command.env
            #
            echo FOO="${FOO[@]}" >> .command.env
            echo /FOO/ >> .command.env
            #
            nxf_eval_cmd STDOUT STDERR bash -c "this --cmd"
            status=$?
            if [ $status -eq 0 ]; then
              echo THIS="$STDOUT" >> .command.env
              echo /THIS/=exit:0 >> .command.env
            else
              echo THIS="$STDERR" >> .command.env
              echo /THIS/=exit:$status >> .command.env
            fi
            #
            nxf_eval_cmd STDOUT STDERR bash -c "other \\"quoted\\" --cmd"
            status=$?
            if [ $status -eq 0 ]; then
              echo THAT="$STDOUT" >> .command.env
              echo /THAT/=exit:0 >> .command.env
            else
              echo THAT="$STDERR" >> .command.env
              echo /THAT/=exit:$status >> .command.env
            fi
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

    def 'should get unstage control script'(){
        given:
        BashWrapperBuilder builder
        when:
        builder = newBashWrapperBuilder()
        then:
        builder.getUnstageControls() == '''\
                cp .command.out /work/dir/.command.out || true
                cp .command.err /work/dir/.command.err || true
                '''.stripIndent()


        when:
        builder = newBashWrapperBuilder(statsEnabled: true)
        then:
        builder.getUnstageControls() == '''\
                cp .command.out /work/dir/.command.out || true
                cp .command.err /work/dir/.command.err || true
                cp .command.trace /work/dir/.command.trace || true
                '''.stripIndent()

        when:
        builder = newBashWrapperBuilder(outputEnvNames: ['some-data'])
        then:
        builder.getUnstageControls() == '''\
                cp .command.out /work/dir/.command.out || true
                cp .command.err /work/dir/.command.err || true
                cp .command.env /work/dir/.command.env || true
                '''.stripIndent()

        when:
        builder = newBashWrapperBuilder(outputEvals: [some:'data'])
        then:
        builder.getUnstageControls() == '''\
                cp .command.out /work/dir/.command.out || true
                cp .command.err /work/dir/.command.err || true
                cp .command.env /work/dir/.command.env || true
                '''.stripIndent()
    }


    def 'should create wrapper with podman' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerEnabled: true,
                containerConfig: [engine: 'podman', enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'podman run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.cleanup_cmd == 'podman rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'podman stop $NXF_BOXID'
    }

    def 'should create wrapper with podman with legacy entrypoint' () {
        when:
        def binding = newBashWrapperBuilder(
                containerImage: 'busybox',
                containerEnabled: true,
                containerConfig: [engine: 'podman', enabled: true, entrypointOverride: true] ).makeBinding()

        then:
        binding.launch_cmd == 'podman run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -w "$NXF_TASK_WORKDIR" --entrypoint /bin/bash --name $NXF_BOXID busybox -c "/bin/bash -ue /work/dir/.command.sh"'
        binding.cleanup_cmd == 'podman rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'podman stop $NXF_BOXID'
    }

    def 'should create wrapper with podman and scratch' () {
        when:
        def binding = newBashWrapperBuilder(
                scratch: true,
                containerImage: 'busybox',
                containerEnabled: true,
                containerConfig: [engine: 'podman', enabled: true] ).makeBinding()

        then:
        binding.launch_cmd == 'podman run -i -e "NXF_TASK_WORKDIR" -v /work/dir:/work/dir -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID busybox /bin/bash -ue /work/dir/.command.sh'
        binding.cleanup_cmd == 'rm -rf $NXF_SCRATCH || true\npodman rm $NXF_BOXID &>/dev/null || true\n'
        binding.kill_cmd == 'podman stop $NXF_BOXID'
    }

    @Unroll
    def 'should check retryable errors' () {
        expect:
        BashWrapperBuilder.isRetryable0(ERROR) == EXPECTED
        where:
        ERROR                           | EXPECTED
        new RuntimeException()          | true
        new SocketException()           | true
        new FileSystemException('foo')  | true
        new IOException()               | false
        new Exception()                 | false
    }
}
