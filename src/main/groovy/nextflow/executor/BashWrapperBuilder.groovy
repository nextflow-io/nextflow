/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.container.ContainerBuilder
import nextflow.container.ContainerScriptTokens
import nextflow.container.DockerBuilder
import nextflow.container.ShifterBuilder
import nextflow.container.SingularityBuilder
import nextflow.container.UdockerBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Escape

/**
 * Builder to create the BASH script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BashWrapperBuilder {

    static final private ENDL = '\n'

    static final public List<String> BASH

    static private level = 0

    @PackageScope
    static public String systemOsName = System.getProperty('os.name')

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
        BASH = Collections.unmodifiableList(  level > 0 ? ['/bin/bash','-uex'] : ['/bin/bash','-ue'] )

    }

    /*
     * Read more linux process handling and signal trap
     *
     * http://mywiki.wooledge.org/SignalTrap
     * http://mywiki.wooledge.org/ProcessManagement
     */
    @PackageScope
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

        nxf_mktemp() {
            local base=${1:-/tmp}
            if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
            else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
            fi
        }
        
        on_exit() {
          exit_status=${ret:=$?}
          printf $exit_status __EXIT_FILE__
          set +u
          [[ "$tee1" ]] && kill $tee1 2>/dev/null
          [[ "$tee2" ]] && kill $tee2 2>/dev/null
          [[ "$ctmp" ]] && rm -rf $ctmp || true
          __EXIT_CMD__
        }

        on_term() {
            set +e
            __KILL_CMD__
        }

        trap on_exit EXIT
        trap on_term TERM INT USR1 USR2
        '''.stripIndent().leftTrim()

    @PackageScope
    static final String SCRIPT_LEGACY = '''
        nxf_tree() {
            declare -a ALL_CHILD
            while read P PP;do
                ALL_CHILD[$PP]+=" $P"
            done < <(ps -e -o pid= -o ppid=)

            stat() {
                local x_ps=$(ps -o pid= -o state= -o pcpu= -o pmem= -o vsz= -o rss= $1)
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
              printf "$tot\\n" || true
            fi
        }

        nxf_sleep() {
          if [[ $1 < 0 ]]; then sleep 5;
          elif [[ $1 < 10 ]]; then sleep 0.1 2>/dev/null || sleep 1;
          elif [[ $1 < 130 ]]; then sleep 1;
          else sleep 5; fi
        }

        nxf_date() {
            local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
        }

        nxf_trace() {
          local pid=$1; local trg=$2;
          local tot;
          local count=0;
          declare -a max=(); for i in {0..13}; do max[i]=0; done
          while [[ true ]]; do
            if ! kill -0 $pid 2>/dev/null; then exit 0; fi
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

    // see http://stackoverflow.com/a/4421282/395921
    //     https://github.com/nextflow-io/nextflow/issues/499
    // resources record:
    // 1    2      3     4     5      5    6       7      8      9      10     11     12          13           14
    // PID, STATE, %CPU, %MEM, VSIZE, RSS, VmPeak, VmHWM, RCHAR, WCHAR, syscr, syscw, READ_BYTES, WRITE_BYTES, CANCEL_W_BYTES
    @PackageScope
    static final String SCRIPT_TRACE = '''   
        set -o pipefail
        prev_total=0
        declare -a prev_time
        mem_tot=$(< /proc/meminfo grep MemTotal | awk '{print $2}')
        num_cpus=$(< /proc/cpuinfo grep '^processor' -c)
        
        nxf_pcpu() {
            local pid=$1
            local proc_time=$(2> /dev/null < /proc/$pid/stat awk '{sum=$14+$15; printf "%.0f",sum}' || echo -n 0)
            local cpu_usage=$(echo -n $proc_time ${prev_time[pid]:-0} $total_time $prev_total $num_cpus | awk '{ pct=($1-$2)/($3-$4)*$5 *100; printf "%.1f", pct }' )
            prev_time[pid]=$proc_time
            nxf_pcpu_ret=$cpu_usage
        }
        
        nxf_tree() {
            declare -a ALL_CHILD
            while read P PP;do
                ALL_CHILD[$PP]+=" $P"
            done < <(ps -e -o pid= -o ppid=)
        
            stat() {
                nxf_pcpu $1 
                local x_pid=$1
                local x_stat=$(2> /dev/null < /proc/$1/stat awk '{print $3}' || echo -n X)
                local x_pcpu=$nxf_pcpu_ret
                
                local x_vsz=$(2> /dev/null < /proc/$1/stat awk '{printf "%.0f", $23/1024}' || echo -n 0)
                local x_rss=$(2> /dev/null < /proc/$1/status grep VmRSS | awk '{print $2}' || echo -n 0)
                local x_pmem=$(echo $x_rss | awk -v mem_tot=$mem_tot '{printf "%.1f", $1/mem_tot*100}')
        
                local x_io=$(2> /dev/null < /proc/$1/io sed 's/^.*:\\s*//' | tr '\\n' ' ' || echo -n 0)
                local x_vm=$(2> /dev/null < /proc/$1/status egrep 'VmPeak|VmHWM' | sed 's/^.*:\\s*//' | sed 's/[\\sa-zA-Z]*$//' | tr '\\n' ' ' || echo -n 0)
        
                stat_ret+="$x_pid $x_stat $x_pcpu $x_pmem $x_vsz $x_rss"
                if [[ $x_vm ]]; then stat_ret+=" $x_vm"; else stat_ret+=" 0 0"; fi
                if [[ $x_io ]]; then stat_ret+=" $x_io"; fi
                stat_ret+='\\n\'
            }
        
            walk() {
                stat $1 
                for i in ${ALL_CHILD[$1]:=}; do walk $i; done
            }
        
            stat_ret=\'\'
            total_time=$(grep '^cpu ' /proc/stat |awk '{sum=$2+$3+$4+$5+$6+$7+$8+$9+$10; printf "%.0f",sum}')
            walk $1
            prev_total=$total_time
            nxf_tree_ret=$stat_ret  
        }
        
        nxf_pstat() {
            nxf_tree $1
            if [[ "$nxf_tree_ret" ]]; then
              nxf_pstat_ret=$(printf "$nxf_tree_ret" | awk '{ t3+=($3*10); t4+=($4*10); t5+=$5; t6+=$6; t7+=$7; t8+=$8; t9+=$9; t10+=$10; t11+=$11; t12+=$12; t13+=$13; t14+=$14 } END { printf "%d 0 %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\\n", NR,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }')
            else
              nxf_pstat_ret=''  
            fi
        }
        
        nxf_sleep() {
          if [[ $1 < 0 ]]; then sleep 5;
          elif [[ $1 < 10 ]]; then sleep 0.1 2>/dev/null || sleep 1;
          elif [[ $1 < 130 ]]; then sleep 1;
          else sleep 5; fi
        }
        
        nxf_date() {
            local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
        }
        
        nxf_trace() {
          local pid=$1; local trg=$2;
          local count=0;
          declare -a max=(); for i in {0..13}; do max[i]=0; done
          while [[ -d /proc/$pid ]]; do
            nxf_pstat $pid
            if [[ "$nxf_pstat_ret" ]]; then
            IFS=' ' read -a val <<< "$nxf_pstat_ret"; unset IFS
            for i in {0..13}; do
              [ ${val[i]} -gt ${max[i]} ] && max[i]=${val[i]}
            done
            echo "pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes" > $trg
            echo "${max[@]}" >> $trg
            fi
            nxf_sleep $count
            count=$((count+1))
          done
        } 

        '''.stripIndent().leftTrim()

    @PackageScope
    static final String MODULE_LOAD = '''
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

        '''.stripIndent().trim()

    @Delegate
    ScriptFileCopyStrategy copyStrategy

    @Delegate
    private TaskBean bean

    private boolean runWithContainer

    private ContainerBuilder containerBuilder

    private Path scriptFile

    private Path inputFile

    private Path startedFile

    private Path exitedFile

    private Path wrapperFile

    private Path stubFile

    BashWrapperBuilder( TaskRun task ) {
        this(new TaskBean(task))
    }

    BashWrapperBuilder( TaskBean bean, ScriptFileCopyStrategy strategy = null ) {
        this.bean = bean
        this.copyStrategy = strategy ?: new SimpleFileCopyStrategy(bean)
    }

    /** only for testing -- do not use */
    protected BashWrapperBuilder() { }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String getScratchDirectoryCommand() {

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
            return 'NXF_SCRATCH="$(set +u; nxf_mktemp $TMPDIR)"'
        }

        if( scratchStr.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NXF_SCRATCH="$(nxf_mktemp /dev/shm)"'
        }

        return "NXF_SCRATCH=\"\$(set +u; nxf_mktemp $scratchStr)\""
    }

    protected boolean fixOwnership() {
        systemOsName == 'Linux' && containerConfig?.fixOwnership && runWithContainer && containerConfig.engine == 'docker' // <-- note: only for docker (shifter is not affected)
    }

    protected Map<String,Path> getResolvedInputs() {
        copyStrategy.resolveForeignFiles(inputFiles)
    }

    /**
     * Build up the BASH wrapper script file which will launch the user provided script
     * @return The {@code Path} of the created wrapper script
     */
    Path build() {
        assert workDir, "Missing 'workDir' property in BashWrapperBuilder object"
        assert script, "Missing 'script' property in BashWrapperBuilder object"

        /* 
         * initialise command files 
         */
        scriptFile = workDir.resolve(TaskRun.CMD_SCRIPT)
        inputFile = workDir.resolve(TaskRun.CMD_INFILE)
        startedFile = workDir.resolve(TaskRun.CMD_START)
        exitedFile = workDir.resolve(TaskRun.CMD_EXIT)
        wrapperFile = workDir.resolve(TaskRun.CMD_RUN)
        stubFile = workDir.resolve(TaskRun.CMD_STUB)
        
        // set true when running with through a container engine
        runWithContainer = containerEnabled && !containerNative

        /*
         * save the input when required
         */
        if( input != null ) {
            Files.write(inputFile, input.toString().getBytes())
        }

        // whenever it has to change to the scratch directory
        final changeDir = getScratchDirectoryCommand()

        /*
         * process the task script
         */
        ContainerScriptTokens scriptTokens = null
        def taskScript = TaskProcessor.normalizeScript(script, shell)
        if(this.containerExecutable) {
            log.warn1("Process `$name` uses an executable container -- This feature has been deprecated and will be removed in a future release", firstOnly: true)
            scriptTokens = ContainerScriptTokens.parse(taskScript)
            environment.putAll( scriptTokens.variables )
        }

        /*
         * create the container launcher command if needed
         */
        containerBuilder = runWithContainer ? createContainerBuilder(changeDir) : null

        /*
         * reformat the task script to include the container launch command when it's a executable container eg:
          *
         * `docker/image:tag --do --something`  ==>  `docker run docker/image:tag --do --something`
         */
        if( scriptTokens ) {
            taskScript = containerBuilder.addContainerRunCommand(scriptTokens)
        }

        /*
         * save the script file
         */
        Files.write(scriptFile, taskScript.getBytes())

        /*
         * fetch the script interpreter i.e. BASH, Perl, Python, etc
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

        wrapper << '# NEXTFLOW TASK: ' << name << ENDL
        wrapper << 'set -e' << ENDL
        wrapper << 'set -u' << ENDL
        wrapper << 'NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x' << ENDL << ENDL

        if( runWithContainer ) {
            containerBuilder.appendHelpers(wrapper)
            // append the process - or - container kill command
            wrapper << scriptCleanUp(exitedFile, changeDir) << ENDL
            // create a random string id to be used an container name
            wrapper << 'export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"' << ENDL
        }
        else {
            wrapper << scriptCleanUp(exitedFile, changeDir) << ENDL
        }

        wrapper << ( changeDir ?: "NXF_SCRATCH=''" ) << ENDL

        // -- print the current environment when debug is enabled
        wrapper << '[[ $NXF_DEBUG > 0 ]] && nxf_env' << ENDL

        def beforeStart = copyStrategy.getBeforeStartScript()
        if( beforeStart ) {
            wrapper << beforeStart << ENDL
        }

        // -- start creating a file to signal that task has began
        wrapper << touchFile(startedFile) << ENDL

        def customEnv = createCustomEnvironment()
        if( customEnv ) {
            // -- enable unbound variables
            wrapper << 'set +u' << ENDL
            wrapper << customEnv 
            // -- disable unbound variables
            wrapper << 'set -u' << ENDL
        }

        /*
         * add the task environment
         */
        def envSnippet = createTaskEnvironment(environment)
        if( envSnippet )
            wrapper <<  envSnippet << ENDL

        /*
         * change in the scratch directory
         */
        wrapper << '[[ $NXF_SCRATCH ]] && echo "nxf-scratch-dir $HOSTNAME:$NXF_SCRATCH" && cd $NXF_SCRATCH' << ENDL

        /*
         * staging input files when required
         */
        def stagingScript = copyStrategy.getStageInputFilesScript(resolvedInputs)
        if( stagingScript ) {
            wrapper << '# stage input files' << ENDL
            wrapper << stagingScript << ENDL
        }

        /*
         * execute the command script
         */
        wrapper << '' << ENDL
        wrapper << 'set +e' << ENDL  // <-- note: use loose error checking so that ops after the script command are executed in all cases
        wrapper << 'ctmp=$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp $TMPDIR)' << ENDL
        wrapper << 'cout=$ctmp/.command.out; mkfifo $cout' << ENDL
        wrapper << 'cerr=$ctmp/.command.err; mkfifo $cerr' << ENDL
        wrapper << 'tee '<< TaskRun.CMD_OUTFILE <<' < $cout &' << ENDL
        wrapper << 'tee1=$!' << ENDL
        wrapper << 'tee '<< TaskRun.CMD_ERRFILE <<' < $cerr >&2 &' << ENDL
        wrapper << 'tee2=$!' << ENDL
        wrapper << '(' << ENDL

        wrapper << getLauncherScript(interpreter,envSnippet) << ENDL

        wrapper << ') >$cout 2>$cerr &' << ENDL
        wrapper << 'pid=$!' << ENDL
        wrapper << 'wait $pid || ret=$?' << ENDL
        wrapper << 'wait $tee1 $tee2' << ENDL

        /*
         * un-stage output files
         */
        if( changeDir ) {
            wrapper << copyFileToWorkDir(TaskRun.CMD_OUTFILE) << ' || true' << ENDL
            wrapper << copyFileToWorkDir(TaskRun.CMD_ERRFILE) << ' || true' << ENDL
        }

        def unstagingScript
        if( (changeDir || workDir != targetDir) && (unstagingScript=copyStrategy.getUnstageOutputFilesScript(outputFiles,targetDir)) ) {
            wrapper << '# copies output files to target' << ENDL
            wrapper << 'if [[ ${ret:=0} == 0 ]]; then' << ENDL
            wrapper << unstagingScript.trim().indent('  ')
            wrapper << 'fi' << ENDL
        }

        if( changeDir && statsEnabled )
            wrapper << copyFileToWorkDir(TaskRun.CMD_TRACE) << ' || true' << ENDL

        if( afterScript ) {
            wrapper << '# user `afterScript`' << ENDL
            wrapper << afterScript << ENDL
        }

        wrapperScript = wrapper.toString()
        Files.write(wrapperFile, wrapperScript.getBytes())
        return wrapperFile
    }

    protected String createCustomEnvironment() {

        def result = new StringBuilder()
        if( beforeScript ) {
            result << '# user `beforeScript`' << ENDL
            result << beforeScript << ENDL
        }

        /*
         * add modules to the environment file
         * note: singularity engine can be loaded by using module
         */
        if( moduleNames ) {
            result << MODULE_LOAD << ENDL << ENDL
            moduleNames.each {
                result << moduleLoad(it) << ENDL
            }
            result << ENDL
        }

        /*
         * activate conda environment
         */
        if( condaEnv ) {
            result << '# conda environment' << ENDL
            result << 'source activate ' << Escape.path(condaEnv) << ENDL
        }

        return result.toString()
    }

    protected String createTaskEnvironment(Map environment) {
        def result = ''
        // when the task needs to be executed by a container
        // wrap the environment in a variable that is passed to the container run
        def wrapper = runWithContainer ? 'nxf_taskenv' : null
        def envSnippet = copyStrategy.getEnvScript(environment, wrapper)
        if( !envSnippet )
            return null

        result += '# task environment' + ENDL
        result += envSnippet
        return result
    }

    protected String getLauncherScript(String interpreter, String env) {

        /*
         * process stats
         */
        String launcher
        final useStubFile = statsEnabled || fixOwnership()
        if( useStubFile ) {
            // create the launcher stub
            createStubScript(interpreter)
            // executes the stub which in turn executes the target command
            launcher = "/bin/bash ${fileStr(stubFile)}"
        }
        else {
            launcher = "${interpreter} ${fileStr(scriptFile)}"
        }

        /*
         * create the container engine command when needed
         */
        if( containerBuilder && !this.containerExecutable) {
            def cmd = env ? 'eval $(nxf_taskenv); ' + launcher : launcher
            launcher = containerBuilder.getRunCommand(cmd)
        }

        /*
         * pipe the input file on the command standard input
         */
        if( !useStubFile && input != null ) {
            launcher += pipeInputFile(inputFile)
        }

        return launcher
    }

    private String copyFileToWorkDir(String fileName) {
        copyFile(fileName, workDir.resolve(fileName))
    }

    private void newLine( StringBuilder builder, String line ) {
        builder.append(line).append(ENDL)
    }

    protected isMacOS() {
        systemOsName.startsWith('Mac')
    }

    protected boolean isLegacyStubScript() {
        isMacOS() && !isContainerEnabled()
    }

    protected Path createStubScript(String interpreter) {
        final stub = new StringBuilder()
        newLine stub, '#!/bin/bash'
        newLine stub, 'set -e'
        newLine stub, 'set -u'
        newLine stub, 'NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 2 ]] && set -x'
        newLine stub, ''
        newLine stub, ( isLegacyStubScript() ? SCRIPT_LEGACY : SCRIPT_TRACE )
        newLine stub, 'trap \'exit ${ret:=$?}\' EXIT'
        newLine stub, 'touch ' + TaskRun.CMD_TRACE
        newLine stub, 'start_millis=$(nxf_date)'
        newLine stub, '('
        String command = "${interpreter} ${fileStr(scriptFile)}"
        if( input != null )
            command += pipeInputFile(inputFile)
        newLine stub, command
        newLine stub, ') &'
        newLine stub, 'pid=$!'                      // get the PID of the main job
        newLine stub, 'nxf_trace "$pid" ' + TaskRun.CMD_TRACE + ' &'
        newLine stub, 'mon=$!'                      // get the pid of the monitor process
        newLine stub, 'wait $pid || ret=$?'         // wait for main job completion and get its exit status
        newLine stub, 'end_millis=$(nxf_date)'      // get the ending time
        newLine stub, 'wait $mon'                   // wait for stats monitor termination
        newLine stub, 'echo $((end_millis-start_millis)) >> ' + TaskRun.CMD_TRACE

        // append to script file chown command to change `root` owner of files created by docker
        if( this.fixOwnership() )
            newLine stub, "\n# patch root ownership problem of files created with docker\n[ \${NXF_OWNER:=''} ] && chown -fR --from root \$NXF_OWNER ${workDir}/{*,.*} || true"

        // save to file
        Files.write(stubFile, stub.toString().getBytes())
        return stubFile
    }

    /**
     * Define the task clean-up snippet
     *
     * @param file The file where the exit status is saved
     * @param dockerKill The command string to kill a container when the task is executed through Docker
     * @return The script string to be included the in main launcher script
     */
    @PackageScope
    String scriptCleanUp( Path file, String scratch ) {

        final script = []

        // -- the kill command
        final killCommand = containerBuilder?.getKillCommand() ?: '[[ "$pid" ]] && nxf_kill $pid'
        // -- cleanup the scratch dir
        if( scratch && cleanup != false ) {
            script << (!containerBuilder ? 'rm -rf $NXF_SCRATCH || true' : '(sudo -n true && sudo rm -rf "$NXF_SCRATCH" || rm -rf "$NXF_SCRATCH")&>/dev/null || true')
        }
        // -- remove the container in this way because 'docker run --rm'  fail in some cases -- see https://groups.google.com/d/msg/docker-user/0Ayim0wv2Ls/-mZ-ymGwg8EJ
        final remove = containerBuilder?.getRemoveCommand()
        if( remove ) {
            script << "${remove} &>/dev/null || true"
        }
        // -- return the exit code
        script << 'exit $exit_status'

        // -- finally compose the script
        SCRIPT_CLEANUP
                .replace('__EXIT_FILE__', exitFile(file))
                .replace('__KILL_CMD__', killCommand)
                .replace('__EXIT_CMD__', script.join('\n  '))
    }


    /**
     * Build a {@link DockerBuilder} object to handle Docker commands
     *
     * @param envFile A file containing environment configuration
     * @param changeDir String command to change to the working directory
     * @return A {@link DockerBuilder} instance
     */
    @PackageScope
    ContainerBuilder createContainerBuilder(String changeDir) {

        final engine = containerConfig.getEngine()
        ContainerBuilder builder

        /*
         * create a builder instance given the container engine
         */
        if( engine == 'docker' )
            builder = new DockerBuilder(containerImage)
        else if( engine == 'singularity' )
            builder = new SingularityBuilder(containerImage)
        else if( engine == 'udocker' )
            builder = new UdockerBuilder(containerImage)
        else if( engine == 'shifter' )
            builder = new ShifterBuilder(containerImage)
        else
            throw new IllegalArgumentException("Unknown container engine: $engine")

        /*
         * initialise the builder
         */
        builder.addMountForInputs(resolvedInputs)

        if( !this.containerExecutable )
            builder.addMount(binDir)

        if(this.containerMount)
            builder.addMount(containerMount)

        // task work dir
        builder.setWorkDir(workDir)

        // set the name
        builder.setName('$NXF_BOXID')

        if( this.containerMemory )
            builder.setMemory(containerMemory)

        if( this.containerCpuset )
            builder.addRunOptions(containerCpuset)

        // export the nextflow script debug variable
        if( statsEnabled || fixOwnership() )
            builder.addEnv( 'NXF_DEBUG=${NXF_DEBUG:=0}')

        // add the user owner variable in order to patch root owned files problem
        if( fixOwnership() )
            builder.addEnv( 'NXF_OWNER=$(id -u):$(id -g)' )

        // set the environment
        if( this.environment && this.containerExecutable ) {
            // PATH variable cannot be extended in an executable container
            // make sure to not include it to avoid to override the container PATH
            environment.remove('PATH')

            // NOTE: the task environment is added only for executable container
            // for *plain* container the environment is passed through the special
            // `nxf_taskenv` bash function wrapper
            builder.addEnv( environment )
        }

        if( engine=='docker' && System.getenv('NXF_DOCKER_OPTS') ) {
            builder.addRunOptions(System.getenv('NXF_DOCKER_OPTS'))
        }

        for( String var : containerConfig.getEnvWhitelist() ) {
            builder.addEnv(var)
        }

        // set up run docker params
        builder.params(containerConfig)

        // extra rule for the 'auto' temp dir temp dir
        def temp = containerConfig.temp?.toString()
        if( temp == 'auto' || temp == 'true' ) {
            builder.setTemp( changeDir ? '$NXF_SCRATCH' : '$(nxf_mktemp)' )
        }

        if( containerConfig.containsKey('kill') )
            builder.params(kill: containerConfig.kill)

        if( containerConfig.writableInputMounts==false )
            builder.params(readOnlyInputs: true)

        // override the docker entry point the image is NOT defined as executable
        if( !this.containerExecutable )
            builder.params(entry: '/bin/bash')

        // give a chance to override any option with process specific `containerOptions`
        if( containerOptions ) {
            builder.addRunOptions(containerOptions)
        }

        builder.build()
        return builder
    }

    String moduleLoad(String name) {
        int p = name.lastIndexOf('/')
        p != -1 ? "nxf_module_load ${name.substring(0,p)} ${name.substring(p+1)}" : "nxf_module_load ${name}"
    }

}
