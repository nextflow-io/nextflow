/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor
import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.DockerBuilder
import nextflow.util.MemoryUnit

/**
 * Builder to create the BASH script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BashWrapperBuilder {

    static final List<String> BASH = ['/bin/bash','-ue']

    /*
     * Read more linux process handling and signal trap
     *
     * http://mywiki.wooledge.org/SignalTrap
     * http://mywiki.wooledge.org/ProcessManagement
     */
    static final String SCRIPT_CLEANUP = '''
        nxf_kill() {
            declare -a ALL_CHILD
            while read P PP;do
                ALL_CHILD[$PP]+=" $P"
            done < <(ps -e -o pid= -o ppid=)

            walk() {
                [[ $1 != $$ ]] && kill $1 2>/dev/null || true
                for i in ${ALL_CHILD[$1]:=}; do walk $i; done
            }

            walk $1
        }

        on_exit() {
          exit_status=${ret:=$?}
          printf $exit_status > __EXIT_FILE__
          exit $exit_status
        }

        on_term() {
            set +e
            __KILL_CMD__
        }

        trap on_exit EXIT
        trap on_term TERM INT USR1 USR2
        '''.stripIndent().leftTrim()

    // note: keep local 'xps' declaration separated by the variable assignment
    // see http://stackoverflow.com/a/4421282/395921
    static final String SCRIPT_TRACE = '''
        nxf_tree() {
            declare -a ALL_CHILD
            while read P PP;do
                ALL_CHILD[$PP]+=" $P"
            done < <(ps -e -o pid= -o ppid=)

            stat() {
                local x_ps=$(ps -o pid=,state=,pcpu=,pmem=,vsz=,rss= $1)
                local x_io=$(cat /proc/$1/io 2> /dev/null | sed 's/^.*:\\s*//' | tr '\\n' ' ')
                local x_vm=$(cat /proc/$1/status 2> /dev/null | egrep 'VmPeak|VmHWM' | sed 's/^.*:\\s*//' | sed 's/[\\sa-zA-Z]*$//' | tr '\\n' ' ')
                [[ ! $x_ps ]] && return 0

                printf "$x_ps"
                if [[ $x_vm ]]; then printf " $x_vm"; else printf " 0 0"; fi
                if [[ $x_io ]]; then printf " $x_io"; fi
                printf "\\n"
            }

            walk() {
                stat $1
                for i in ${ALL_CHILD[$1]:=}; do walk $i; done
            }

            walk $1
        }

        nxf_pstat() {
            local data=$(nxf_tree $1)
            local tot=\'\'
            if [[ "$data" ]]; then
              tot=$(awk '{ t3+=($3*10); t4+=($4*10); t5+=$5; t6+=$6; t7+=$7; t8+=$8; t9+=$9; t10+=$10; t11+=$11; t12+=$12; t13+=$13; t14+=$14 } END { print NR,"0",t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }' <<< "$data")
              printf "$tot\\n"
            fi
        }

        nxf_sleep() {
          if [[ $1 < 0 ]]; then sleep 5;
          elif [[ $1 < 10 ]]; then sleep 0.1;
          elif [[ $1 < 130 ]]; then sleep 1;
          else sleep 5; fi
        }

        nxf_date() {
            case `uname` in
                Darwin) if hash gdate 2>/dev/null; then echo 'gdate +%s%3N'; else echo 'date +%s000'; fi;;
                *) echo 'date +%s%3N';;
            esac
        }

        NXF_DATE=$(nxf_date)

        nxf_trace() {
          local pid=$1; local trg=$2;
          local tot;
          local count=0;
          declare -a max=(); for i in {0..13}; do max[i]=0; done
          while [[ true ]]; do
            tot=$(nxf_pstat $pid)
            [[ ! $tot ]] && break
            IFS=' ' read -a val <<< "$tot"; unset IFS
            for i in {0..13}; do
              [ ${val[i]} -gt ${max[i]} ] && max[i]=${val[i]}
            done
            echo "pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes" > $trg
            echo "${max[@]}" >> $trg
            nxf_sleep $count
            count=$((count+1))
          done
        }

        '''.stripIndent().leftTrim()

    final private TaskRun task

    def scratch

    def input

    Map<String,String> environment

    String headerScript

    String stagingScript

    String unstagingScript

    String wrapperScript

    String dockerImage

    List<String> moduleNames

    Path workDir

    Path targetDir

    String script

    def shell

    Map dockerConfig

    Path dockerMount

    String dockerCpuset

    MemoryUnit dockerMemory

    String name

    private runWithDocker

    boolean statsEnabled

    String beforeScript

    String afterScript

    BashWrapperBuilder( TaskRun task ) {
        this.task = task
        this.name = "nxf-" + task.hash?.toString()?.substring(0,8)

        // set the input (when available)
        this.input = task.stdin
        this.scratch = task.scratch
        this.workDir = task.workDir
        this.targetDir = task.targetDir

        // set the environment
        // note: create a copy of the process environment to avoid concurrent
        // process executions override each others
        this.environment = new HashMap( task.processor.getProcessEnvironment() )
        this.environment.putAll( task.getInputEnvironment() )

        this.moduleNames = task.config.getModule()
        this.shell = task.config.shell
        this.script = task.script.toString()
        this.beforeScript = task.config.beforeScript
        this.afterScript = task.config.afterScript

        // docker config
        this.dockerImage = task.container
        this.dockerConfig = task.processor?.session?.config?.docker
        this.dockerMemory = task.config.getMemory()

        // stats
        this.statsEnabled = task.processor?.session?.statsEnabled
    }

    BashWrapperBuilder( Map params ) {
        log.trace "Wrapper params: $params"

        task = null
        this.name = params.name
        this.shell = params.shell ?: BASH
        this.script = params.script?.toString()
        this.input = params.input
        this.scratch = params.scratch
        this.workDir = params.workDir
        this.targetDir = params.targetDir
        this.environment = params.environment
        this.headerScript = params.headerScript
        this.moduleNames = params.moduleNames
        this.beforeScript = params.beforeScript
        this.afterScript = params.afterScript

        // docker config
        this.dockerImage = params.container
        this.dockerConfig = params.dockerConfig
        this.dockerMount = params.dockerMount
        this.statsEnabled = params.statsEnabled
    }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String changeToScratchDirectory() {

        if( scratch == null || scratch == false ) {
            return null
        }

        /*
         * when 'scratch' is defined as a bool value
         * try to use the 'TMP' variable, if does not exist fallback to a tmp folder
         */
        if( scratch == true ) {
            return 'NXF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NXF_SCRATCH'
        }

        // convert to string for safety
        final scratchStr = scratch.toString()

        // when it is defined by a variable, just use it
        if( scratchStr.startsWith('$') ) {
            return "NXF_SCRATCH=\${${scratchStr.substring(1)}:-`mktemp -d`} && cd \$NXF_SCRATCH"
        }

        if( scratchStr.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NXF_SCRATCH=$(mktemp -d -p /dev/shm/) && cd $NXF_SCRATCH'
        }


        return "NXF_SCRATCH=\$(mktemp -d -p $scratch) && cd \$NXF_SCRATCH"

    }

    /**
     * Build up the BASH wrapper script file which will launch the user provided script
     * @return The {@code Path} of the created wrapper script
     */

    Path build() {
        assert workDir, "Missing 'workDir' property in BashWrapperBuilder object"
        assert script, "Missing 'script' property in BashWrapperBuilder object"

        final ENDL = '\n'
        final scriptFile = workDir.resolve(TaskRun.CMD_SCRIPT)
        final inputFile = workDir.resolve(TaskRun.CMD_INFILE)
        final environmentFile = workDir.resolve(TaskRun.CMD_ENV)
        final startedFile = workDir.resolve(TaskRun.CMD_START)
        final exitedFile = workDir.resolve(TaskRun.CMD_EXIT)
        final runnerFile = workDir.resolve(TaskRun.CMD_RUN)
        final wrapperFile = workDir.resolve(TaskRun.CMD_WRAPPER)

        // set true when running with docker
        runWithDocker = dockerImage && dockerConfig?.enabled?.toString() == 'true'

        /*
         * the script file
         */
        final taskScript = scriptFile.text = TaskProcessor.normalizeScript(script, shell)

        /*
         * fetch the script interpreter
         */
        final interpreter = TaskProcessor.fetchInterpreter(taskScript)

        /*
         * save the input when required
         */
        if( input != null ) {
            inputFile.text = input
        }

        /*
         * add modules to the environment file
         */
        moduleNames?.each { String name ->
            environmentFile << "module load $name" << ENDL
        }

        /*
         * save the environment
         */
        if( environment ) {
            // create the *bash* environment script
            environmentFile << TaskProcessor.bashEnvironmentScript(environment)
        }

        // whenever it has to change to the scratch directory
        final changeDir = changeToScratchDirectory()

        DockerBuilder docker = runWithDocker ? createDockerBuilder(environmentFile,changeDir) : null

        /*
         * create a script runner which do the following
         * 1 - move the TMP directory provided by the sge/oge grid engine
         * 2 - pipe the input stream
         * 3 - launch the user script
         * 4 - un-stage e.g. copy back the result files to the working folder
         */

        /*
         * Note: The signals SIGKILL and SIGSTOP cannot be caught, blocked, or ignored.
         * Read more: http://man7.org/linux/man-pages/man7/signal.7.html
         */

        def runner = new StringBuilder()
        runner << '#!' << BASH.join(' ') << ENDL
        if( headerScript )
            runner << headerScript << ENDL

        runner << scriptCleanUp(exitedFile, docker?.killCommand) << ENDL
        runner << touchFile(startedFile) << ENDL

        if( beforeScript ) {
            runner << '# user `beforeScript`' << ENDL
            runner << beforeScript << ENDL
        }

        // source the environment
        if( !runWithDocker ) {
            runner << '[ -f '<< fileStr(environmentFile) << ' ]' << ' && source ' << fileStr(environmentFile) << ENDL
        }

        if( changeDir ) {
            runner << changeDir << ENDL
        }

        // staging input files when required
        if( stagingScript ) {
            runner << stagingScript << ENDL
        }

        // execute the command script
        runner << '' << ENDL
        runner << 'set +e' << ENDL  // <-- note: use loose error checking so that ops after the script command are executed in all cases
        runner << '(' << ENDL

        // execute by invoking the command through a Docker container
        if( docker ) {
            runner << docker.runCommand << ' '
        }

        /*
         * process stats
         */
        if( statsEnabled ) {
            final wrapper = new StringBuilder()
            wrapper << '#!' << BASH.join(' ') << ENDL
            wrapper << SCRIPT_TRACE << ENDL
            wrapper << 'trap \'exit ${ret:=$?}\' EXIT' << ENDL
            wrapper << 'start_millis=$($NXF_DATE)'  << ENDL
            wrapper << '(' << ENDL
            wrapper << interpreter << ' ' << fileStr(scriptFile)
            if( input != null ) wrapper << ' < ' << fileStr(inputFile)
            wrapper << ' &> ' << TaskRun.CMD_OUTFILE << ENDL
            wrapper << ') &' << ENDL
            wrapper << 'pid=$!' << ENDL                     // get the PID of the main job
            wrapper << 'nxf_trace "$pid" ' << TaskRun.CMD_TRACE << ' &' << ENDL
            wrapper << 'mon=$!' << ENDL                     // get the pid of the monitor process
            wrapper << 'wait $pid' << ENDL                  // wait for main job completion
            wrapper << 'ret=$?' << ENDL                     // get the main job error code
            wrapper << 'end_millis=$($NXF_DATE)' << ENDL    // get the ending time
            wrapper << 'kill $mon || wait $mon' << ENDL     // kill the monitor and wait for its ending
            wrapper << '[ -f ' << TaskRun.CMD_TRACE << ' ] && echo $((end_millis-start_millis)) >> ' << TaskRun.CMD_TRACE << ENDL
            // save to file
            wrapperFile.text = wrapper.toString()

            // invoke it from the main script
            runner << BASH.join(' ') << ' ' << fileStr(wrapperFile) << ENDL
        }
        else {
            runner << interpreter << ' ' << fileStr(scriptFile)
            if( input != null ) runner << ' < ' << fileStr(inputFile)
            runner << ' &> ' << TaskRun.CMD_OUTFILE << ENDL
        }
        runner << ') &' << ENDL
        runner << 'pid=$!' << ENDL
        runner << 'wait $pid || ret=$?' << ENDL

        /*
         * docker clean-up
         */
        if( docker?.removeCommand ) {
            // remove the container in this way because 'docker run --rm'  fail in some cases -- see https://groups.google.com/d/msg/docker-user/0Ayim0wv2Ls/-mZ-ymGwg8EJ
            runner << docker.removeCommand << ' &>/dev/null &' << ENDL
        }

        /*
         * un-stage output files
         */
        if( changeDir )
            runner << copyFile(TaskRun.CMD_OUTFILE, workDir) << ' || true' << ENDL

        if( (changeDir || workDir != targetDir) && unstagingScript  )
            runner << unstagingScript << ENDL

        if( changeDir && statsEnabled )
            runner << copyFile(TaskRun.CMD_TRACE,  workDir) << ' || true' << ENDL

        if( afterScript ) {
            runner << '# user `afterScript`' << ENDL
            runner << afterScript << ENDL
        }

        runnerFile.text = wrapperScript = runner.toString()
        return runnerFile
    }

    /**
     * Define the task clean-up snippet
     *
     * @param file The file where the exit status is saved
     * @param dockerKill The command string to kill a container when the task is executed through Docker
     * @return The script string to be included the in main launcher script
     */
    @PackageScope
    String scriptCleanUp( Path file, String dockerKill ) {
        SCRIPT_CLEANUP
                .replace('__EXIT_FILE__', exitFile(file))
                .replace('__KILL_CMD__', dockerKill ?: '[[ "$pid" ]] && nxf_kill $pid')
    }

    /**
     * Build a {@link DockerBuilder} object to handle Docker commands
     *
     * @param envFile A file containing environment configuration
     * @param changeDir String command to change to the working directory
     * @return A {@link DockerBuilder} instance
     */
    @PackageScope
    DockerBuilder createDockerBuilder(Path envFile, String changeDir) {

        def docker = new DockerBuilder(dockerImage)
        if( task ) {
            docker.addMountForInputs( task.getInputFiles() )
                    .addMount( task.processor.session.workDir )
                    .addMount( task.processor.session.binDir )
        }

        // set the name
        docker.setName(this.name)

        if( dockerMount )
            docker.addMount(dockerMount)

        if( dockerMemory )
            docker.setMemory(dockerMemory)

        if( dockerCpuset )
            docker.setCpus(dockerCpuset)

        // set the environment
        if( !envFile.empty() )
            docker.addEnv( envFile )

        // turn on container remove by default
        if( !dockerConfig.containsKey('remove') )
            dockerConfig.remove = true

        // set up run docker params
        docker.params(dockerConfig)

        // extra rule for the 'auto' temp dir temp dir
        def temp = dockerConfig.temp?.toString()
        if( temp == 'auto' || temp == 'true' ) {
            docker.setTemp( changeDir ? '$NXF_SCRATCH' : '$(mktemp -d)' )
        }

        docker.build()
        return docker
    }


    protected String touchFile( Path file ) {
        "touch ${file.toString()}"
    }

    protected String copyFile( String name, Path target ) {
        "cp ${name} ${target}"
    }

    protected String exitFile( Path file ) {
        file.toString()
    }

    protected String fileStr( Path file ) {
        file.toString()
    }

}
