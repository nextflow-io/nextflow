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
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.ContainerScriptTokens
import nextflow.util.DockerBuilder
/**
 * Builder to create the BASH script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BashWrapperBuilder {

    static final private ENDL = '\n'

    static final List<String> BASH

    static private level = 0

    @PackageScope
    static String systemOsName = System.getProperty('os.name')

    static {
        /*
         * Env variable `NXF_DEBUG` is used to control debug options in executed BASH scripts
         * - 0: no debug
         * - 1: dump current environment in the `.command.log` file
         * - 2: trace the execution of user script adding the `set -x` flag
         * - 3: trace the execution of wrapper scripts
         */
        def str = System.getenv('NXF_DEBUG')
        try {
            level = str as Integer
        }
        catch( Exception e ) {
            log.warn "Invalid value for `NXF_DEBUG` variable: $str -- See http://www.nextflow.io/docs/latest/config.html#environment-variables"
        }
        BASH = Collections.unmodifiableList(  level > 1 ? ['/bin/bash','-uex'] : ['/bin/bash','-ue'] )

    }

    /*
     * Read more linux process handling and signal trap
     *
     * http://mywiki.wooledge.org/SignalTrap
     * http://mywiki.wooledge.org/ProcessManagement
     */
    static final String SCRIPT_CLEANUP = '''
        nxf_env() {
            echo '============= task environment ============='
            env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
            echo '============= task output =================='
        }

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

        function nxf_mktemp() {
            local base=${1:-/tmp}
            [[ $(uname) = Darwin ]] && mktemp -d $base/nxf.XXXXXXXXXX || mktemp -d -t nxf.XXXXXXXXXX -p $base
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
    // resources record:
    // 1    2      3     4     5      5    6       7      8      9      10     11     12          13           14
    // PID, STATE, %CPU, %MEM, VSIZE, RSS; VmPeak, VmHWM; RCHAR, WCHAR, syscr, syscw, READ_BYTES, WRITE_BYTES, CANCEL_W_BYTES
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
              tot=$(awk '{ t3+=($3*10); t4+=($4*10); t5+=$5; t6+=$6; t7+=$7; t8+=$8; t9+=$9; t10+=$10; t11+=$11; t12+=$12; t13+=$13; t14+=$14 } END { printf "%d 0 %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\\n", NR,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }' <<< "$data")
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

    @Delegate
    ScriptFileCopyStrategy copyStrategy

    @Delegate
    private TaskBean bean

    private runWithDocker

    BashWrapperBuilder( TaskRun task ) {
        this(new TaskBean(task))
    }


    BashWrapperBuilder( TaskBean bean, ScriptFileCopyStrategy strategy = null ) {
        this.bean = bean
        this.copyStrategy = strategy ?: new SimpleFileCopyStrategy(bean)
    }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String changeToScratchDirectory() {

        // convert to string for safety
        final scratchStr = scratch?.toString()

        if( scratchStr == null || scratchStr == 'false' ) {
            return null
        }

        /*
         * when 'scratch' is defined as a bool value
         * try to use the 'TMP' variable, if does not exist fallback to a tmp folder
         */
        if( scratchStr == 'true' ) {
            return 'NXF_SCRATCH="$(set +u; nxf_mktemp $TMPDIR)" && cd $NXF_SCRATCH'
        }

        if( scratchStr.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NXF_SCRATCH="$(nxf_mktemp /dev/shm/)" && cd $NXF_SCRATCH'
        }

        return "NXF_SCRATCH=\"\$(set +u; nxf_mktemp $scratchStr)\" && cd \$NXF_SCRATCH"
    }

    /**
     * Build up the BASH wrapper script file which will launch the user provided script
     * @return The {@code Path} of the created wrapper script
     */

    Path build() {
        assert workDir, "Missing 'workDir' property in BashWrapperBuilder object"
        assert script, "Missing 'script' property in BashWrapperBuilder object"

        final scriptFile = workDir.resolve(TaskRun.CMD_SCRIPT)
        final inputFile = workDir.resolve(TaskRun.CMD_INFILE)
        final environmentFile = workDir.resolve(TaskRun.CMD_ENV)
        final startedFile = workDir.resolve(TaskRun.CMD_START)
        final exitedFile = workDir.resolve(TaskRun.CMD_EXIT)
        final wrapperFile = workDir.resolve(TaskRun.CMD_RUN)
        final stubFile = workDir.resolve(TaskRun.CMD_STUB)

        // set true when running with docker
        runWithDocker = dockerImage && (executable || dockerConfig.enabled?.toString() == 'true')

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

        /*
         * process the task script
         */
        ContainerScriptTokens scriptTokens = null
        def taskScript = TaskProcessor.normalizeScript(script, shell)
        if( executable ) {
            scriptTokens = ContainerScriptTokens.parse(taskScript)
            environment.putAll( scriptTokens.variables )
        }

        /*
         * create the docker command if required
         */
        DockerBuilder docker = runWithDocker ? createDockerBuilder(environment,changeDir) : null

        /*
         * save the script file
         */
        if( scriptTokens ) {
            taskScript = addContainerRunCommand( scriptTokens, docker )
        }
        scriptFile.text = taskScript

        /*
         * fetch the script interpreter
         */
        final interpreter = TaskProcessor.fetchInterpreter(taskScript)

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

        def wrapper = new StringBuilder()
        wrapper << '#!/bin/bash' << ENDL
        if( headerScript )
            wrapper << headerScript << ENDL

        wrapper << 'set -e' << ENDL
        wrapper << 'set -u' << ENDL
        wrapper << 'NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 2 ]] && set -x' << ENDL << ENDL

        if( runWithDocker ) {
            // append the process - or - container kill command
            wrapper << scriptCleanUp(exitedFile, docker.killCommand) << ENDL
            // create a random string id to be used an container name
            wrapper << 'export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"' << ENDL
            // append to script file chown command to change `root` owner of files created by docker
            // note: this is not required when running docker in OSX through boo2docker because it manage correctly files ownership
            if( systemOsName == 'Linux' && dockerConfig.fixOwnership )
            scriptFile << "\n# patch root ownership problem of files created with docker\n[ \${NXF_OWNER:=''} ] && chown -fR --from root \$NXF_OWNER ${workDir}/{*,.*} || true\n"
        }
        else {
            wrapper << scriptCleanUp(exitedFile) << ENDL
        }

        // -- print the current environment when debug is enabled
        wrapper << '[[ $NXF_DEBUG > 0 ]] && nxf_env' << ENDL

        def beforeStart = copyStrategy.getBeforeStartScript()
        if( beforeStart ) {
            wrapper << beforeStart << ENDL
        }

        // -- start creating a file to signal that task has began
        wrapper << touchFile(startedFile) << ENDL

        if( beforeScript ) {
            wrapper << '# user `beforeScript`' << ENDL
            wrapper << beforeScript << ENDL
        }

        // source the environment
        if( !runWithDocker ) {
            wrapper << '[ -f '<< fileStr(environmentFile) << ' ]' << ' && source ' << fileStr(environmentFile) << ENDL
        }

        if( changeDir ) {
            wrapper << changeDir << ENDL
        }

        // staging input files when required
        def stagingScript = copyStrategy.getStageInputFilesScript()
        if( stagingScript ) {
            wrapper << stagingScript << ENDL
        }

        // execute the command script
        wrapper << '' << ENDL
        wrapper << 'set +e' << ENDL  // <-- note: use loose error checking so that ops after the script command are executed in all cases
        wrapper << 'COUT=$PWD/.command.po; mkfifo $COUT' << ENDL
        wrapper << 'CERR=$PWD/.command.pe; mkfifo $CERR' << ENDL
        wrapper << 'tee '<< TaskRun.CMD_OUTFILE <<' < $COUT &' << ENDL
        wrapper << 'tee1=$!' << ENDL
        wrapper << 'tee '<< TaskRun.CMD_ERRFILE <<' < $CERR >&2 &' << ENDL
        wrapper << 'tee2=$!' << ENDL
        wrapper << '(' << ENDL

        // execute by invoking the command through a Docker container
        if( docker && !executable ) {
            wrapper << docker.runCommand << " -c '"
        }

        /*
         * process stats
         */
        if( statsEnabled ) {
            final stub = new StringBuilder()
            stub << '#!/bin/bash' << ENDL
            stub << 'set -e' << ENDL
            stub << 'set -u' << ENDL
            stub << 'NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 3 ]] && set -x' << ENDL << ENDL
            stub << SCRIPT_TRACE << ENDL
            stub << 'trap \'exit ${ret:=$?}\' EXIT' << ENDL
            stub << 'touch ' << TaskRun.CMD_TRACE << ENDL
            stub << 'start_millis=$($NXF_DATE)'  << ENDL
            stub << '(' << ENDL
            stub << interpreter << ' ' << fileStr(scriptFile)
            if( input != null ) stub << ' < ' << fileStr(inputFile)
            stub << ENDL
            stub << ') &' << ENDL
            stub << 'pid=$!' << ENDL                     // get the PID of the main job
            stub << 'nxf_trace "$pid" ' << TaskRun.CMD_TRACE << ' &' << ENDL
            stub << 'mon=$!' << ENDL                     // get the pid of the monitor process
            stub << 'wait $pid || ret=$?' << ENDL        // wait for main job completion and get its exit status
            stub << 'end_millis=$($NXF_DATE)' << ENDL    // get the ending time
            stub << 'kill $mon || wait $mon' << ENDL     // kill the monitor and wait for its ending
            stub << 'echo $((end_millis-start_millis)) >> ' << TaskRun.CMD_TRACE << ENDL
            // save to file
            stubFile.text = stub.toString()

            // invoke it from the main script
            wrapper << '/bin/bash ' << fileStr(stubFile)
            if( docker && !executable ) wrapper << "'"
        }
        else {
            wrapper << interpreter << ' ' << fileStr(scriptFile)
            if( docker && !executable ) wrapper << "'"
            if( input != null ) wrapper << ' < ' << fileStr(inputFile)
        }
        wrapper << ENDL
        wrapper << ') >$COUT 2>$CERR &' << ENDL
        wrapper << 'pid=$!' << ENDL
        wrapper << 'wait $pid || ret=$?' << ENDL
        wrapper << 'wait $tee1 $tee2' << ENDL

        /*
         * docker clean-up
         */
        if( docker?.removeCommand ) {
            // remove the container in this way because 'docker run --rm'  fail in some cases -- see https://groups.google.com/d/msg/docker-user/0Ayim0wv2Ls/-mZ-ymGwg8EJ
            wrapper << docker.removeCommand << ' &>/dev/null &' << ENDL
        }

        /*
         * un-stage output files
         */
        if( changeDir ) {
            wrapper << copyFile(TaskRun.CMD_OUTFILE, workDir) << ' || true' << ENDL
            wrapper << copyFile(TaskRun.CMD_ERRFILE, workDir) << ' || true' << ENDL
        }

        def unstagingScript
        if( (changeDir || workDir != targetDir) && (unstagingScript=copyStrategy.getUnstageOutputFilesScript()) )
            wrapper << unstagingScript << ENDL

        if( changeDir && statsEnabled )
            wrapper << copyFile(TaskRun.CMD_TRACE,  workDir) << ' || true' << ENDL

        if( afterScript ) {
            wrapper << '# user `afterScript`' << ENDL
            wrapper << afterScript << ENDL
        }

        wrapperFile.text = wrapperScript = wrapper.toString()
        return wrapperFile
    }

    /**
     * Define the task clean-up snippet
     *
     * @param file The file where the exit status is saved
     * @param dockerKill The command string to kill a container when the task is executed through Docker
     * @return The script string to be included the in main launcher script
     */
    @PackageScope
    String scriptCleanUp( Path file, String dockerKill = null ) {
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
    DockerBuilder createDockerBuilder(Map environment, String changeDir) {

        def docker = new DockerBuilder(dockerImage)
        docker.addMountForInputs(inputFiles)
        docker.addMount(workDir)
        if( !executable )
            docker.addMount(binDir)

        // set the name
        docker.setName('$NXF_BOXID')

        if( dockerMemory )
            docker.setMemory(dockerMemory)

        if( dockerCpuset )
            docker.addRunOptions(dockerCpuset)

        // set the environment
        if( environment ) {
            // export the nextflow script debug variable
            docker.addEnv( 'NXF_DEBUG=${NXF_DEBUG:=0}')

            // add the user owner variable in order to patch root owned files problem
            if( dockerConfig.fixOwnership )
                docker.addEnv( 'NXF_OWNER=$(id -u):$(id -g)' )

            if( executable ) {
                // PATH variable cannot be extended in an executable container
                // make sure to not include it to avoid to override the container PATH
                environment.remove('PATH')
                docker.addEnv( environment )
            }
            else
                docker.addEnv( workDir.resolve(TaskRun.CMD_ENV) )

        }

        // set up run docker params
        docker.params(dockerConfig)

        // extra rule for the 'auto' temp dir temp dir
        def temp = dockerConfig.temp?.toString()
        if( temp == 'auto' || temp == 'true' ) {
            docker.setTemp( changeDir ? '$NXF_SCRATCH' : '$(nxf_mktemp)' )
        }

        if( dockerConfig.containsKey('kill') )
            docker.params(kill: dockerConfig.kill)

        // override the docker entry point the image is NOT defined as executable
        if( !executable )
            docker.params(entry: '/bin/bash')

        docker.build()
        return docker
    }



    /**
     * Given a normalised shell script (starting with a she-bang line)
     * replace the first token on the first line with a docker run command
     *
     * @param script
     * @param docker
     * @return
     */
    String addContainerRunCommand( ContainerScriptTokens script, DockerBuilder docker ) {

        final result = new ArrayList<String>(script.lines)
        final i = script.index
        final main = result[i].trim()
        final p = main.indexOf(' ')
        result[i] = ( p != -1
                    ? docker.runCommand + main.substring(p)
                    : docker.runCommand )
        result.add('')
        return result.join('\n')
    }
}
