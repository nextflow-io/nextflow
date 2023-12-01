/*
 * copyright 2013-2023, seqera labs
 *
 * licensed under the apache license, version 2.0 (the "license");
 * you may not use this file except in compliance with the license.
 * you may obtain a copy of the license at
 *
 *     http://www.apache.org/licenses/license-2.0
 *
 * unless required by applicable law or agreed to in writing, software
 * distributed under the license is distributed on an "as is" basis,
 * without warranties or conditions of any kind, either express or implied.
 * see the license for the specific language governing permissions and
 * limitations under the license.
 */

package nextflow.cli

import java.nio.file.nosuchfileexception
import java.nio.file.path
import java.util.regex.pattern

import com.beust.jcommander.dynamicparameter
import com.beust.jcommander.istringconverter
import com.beust.jcommander.parameter
import com.beust.jcommander.parameters
import groovy.json.jsonslurper
import groovy.transform.compilestatic
import groovy.transform.memoized
import groovy.util.logging.slf4j
import groovyx.gpars.gparsconfig
import nextflow.const
import nextflow.nf
import nextflow.nextflowmeta
import nextflow.sysenv
import nextflow.config.configbuilder
import nextflow.config.configmap
import nextflow.exception.abortoperationexception
import nextflow.file.filehelper
import nextflow.plugin.plugins
import nextflow.scm.assetmanager
import nextflow.script.scriptfile
import nextflow.script.scriptrunner
import nextflow.secret.secretsloader
import nextflow.util.custompoolfactory
import nextflow.util.duration
import nextflow.util.historyfile
import org.yaml.snakeyaml.yaml
/**
 * cli sub-command run
 *
 * @author paolo di tommaso <paolo.ditommaso@gmail.com>
 */
@slf4j
@compilestatic
@parameters(commanddescription = "execute a pipeline project")
class cmdrun extends cmdbase implements huboptions {

    static final public pattern run_name_pattern = pattern.compile(/^[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}$/, pattern.case_insensitive)

    static final public list<string> valid_params_file = ['json', 'yml', 'yaml']

    static final public dsl2 = '2'
    static final public dsl1 = '1'

    static {
        // install the custom pool factory for gpars threads
        gparsconfig.poolfactory = new custompoolfactory()
    }

    static class durationconverter implements istringconverter<long> {
        @override
        long convert(string value) {
            if( !value ) throw new illegalargumentexception()
            if( value.islong() ) {  return value.tolong() }
            return duration.of(value).tomillis()
        }
    }

    static final public string name = 'run'

    private map<string,string> sysenv = system.getenv()

    @parameter(names=['-name'], description = 'assign a mnemonic name to the a pipeline run')
    string runname

    @parameter(names=['-lib'], description = 'library extension path')
    string libpath

    @parameter(names=['-cache'], description = 'enable/disable processes caching', arity = 1)
    boolean cacheable

    @parameter(names=['-resume'], description = 'execute the script using the cached results, useful to continue executions that was stopped by an error')
    string resume

    @parameter(names=['-ps','-pool-size'], description = 'number of threads in the execution pool', hidden = true)
    integer poolsize

    @parameter(names=['-pi','-poll-interval'], description = 'executor poll interval (duration string ending with ms|s|m)', converter = durationconverter, hidden = true)
    long pollinterval

    @parameter(names=['-qs','-queue-size'], description = 'max number of processes that can be executed in parallel by each executor')
    integer queuesize

    @parameter(names=['-test'], description = 'test a script function with the name specified')
    string test

    @parameter(names=['-w', '-work-dir'], description = 'directory where intermediate result files are stored')
    string workdir

    @parameter(names=['-bucket-dir'], description = 'remote bucket where intermediate result files are stored')
    string bucketdir

    @parameter(names=['-with-cloudcache'], description = 'enable the use of object storage bucket as storage for cache meta-data')
    string cloudcachepath

    /**
     * defines the parameters to be passed to the pipeline script
     */
    @dynamicparameter(names = '--', description = 'set a parameter used by the pipeline', hidden = true)
    map<string,string> params = new linkedhashmap<>()

    @parameter(names='-params-file', description = 'load script parameters from a json/yaml file')
    string paramsfile

    @dynamicparameter(names = ['-process.'], description = 'set process options' )
    map<string,string> process = [:]

    @dynamicparameter(names = ['-e.'], description = 'add the specified variable to execution environment')
    map<string,string> env = [:]

    @parameter(names = ['-e'], description = 'exports all current system environment')
    boolean exportsysenv

    @dynamicparameter(names = ['-executor.'], description = 'set executor options', hidden = true )
    map<string,string> executoroptions = [:]

    @parameter(description = 'project name or repository url')
    list<string> args

    @parameter(names=['-r','-revision'], description = 'revision of the project to run (either a git branch, tag or commit sha number)')
    string revision

    @parameter(names=['-d','-deep'], description = 'create a shallow clone of the specified depth')
    integer deep

    @parameter(names=['-latest'], description = 'pull latest changes before run')
    boolean latest

    @parameter(names='-stdin', hidden = true)
    boolean stdin

    @parameter(names = ['-ansi'], hidden = true, arity = 0)
    void setansi(boolean value) {
        launcher.options.ansilog = value
    }

    @parameter(names = ['-ansi-log'], description = 'enable/disable ansi console logging', arity = 1)
    void setansilog(boolean value) {
        launcher.options.ansilog = value
    }

    @parameter(names = ['-with-tower'], description = 'monitor workflow execution with seqera tower service')
    string withtower

    @parameter(names = ['-with-wave'], hidden = true)
    string withwave

    @parameter(names = ['-with-fusion'], hidden = true)
    string withfusion

    @parameter(names = ['-with-weblog'], description = 'send workflow status messages via http to target url')
    string withweblog

    @parameter(names = ['-with-trace'], description = 'create processes execution tracing file')
    string withtrace

    @parameter(names = ['-with-report'], description = 'create processes execution html report')
    string withreport

    @parameter(names = ['-with-timeline'], description = 'create processes execution timeline file')
    string withtimeline

    @parameter(names = '-with-charliecloud', description = 'enable process execution in a charliecloud container runtime')
    def withcharliecloud

    @parameter(names = '-with-singularity', description = 'enable process execution in a singularity container')
    def withsingularity

    @parameter(names = '-with-apptainer', description = 'enable process execution in a apptainer container')
    def withapptainer

    @parameter(names = '-with-podman', description = 'enable process execution in a podman container')
    def withpodman

    @parameter(names = '-without-podman', description = 'disable process execution in a podman container')
    def withoutpodman

    @parameter(names = '-with-docker', description = 'enable process execution in a docker container')
    def withdocker

    @parameter(names = '-without-docker', description = 'disable process execution with docker', arity = 0)
    boolean withoutdocker

    @parameter(names = '-with-mpi', hidden = true)
    boolean withmpi

    @parameter(names = '-with-dag', description = 'create pipeline dag file')
    string withdag

    @parameter(names = ['-bg'], arity = 0, hidden = true)
    void setbackground(boolean value) {
        launcher.options.background = value
    }

    @parameter(names=['-c','-config'], hidden = true )
    list<string> runconfig

    @dynamicparameter(names = ['-cluster.'], description = 'set cluster options', hidden = true )
    map<string,string> clusteroptions = [:]

    @parameter(names=['-profile'], description = 'choose a configuration profile')
    string profile

    @parameter(names=['-dump-hashes'], description = 'dump task hash keys for debugging purpose')
    string dumphashes

    @parameter(names=['-dump-channels'], description = 'dump channels for debugging purpose')
    string dumpchannels

    @parameter(names=['-n','-with-notification'], description = 'send a notification email on workflow completion to the specified recipients')
    string withnotification

    @parameter(names=['-with-conda'], description = 'use the specified conda environment package or file (must end with .yml|.yaml suffix)')
    string withconda

    @parameter(names=['-without-conda'], description = 'disable the use of conda environments')
    boolean withoutconda

    @parameter(names=['-with-spack'], description = 'use the specified spack environment package or file (must end with .yaml suffix)')
    string withspack

    @parameter(names=['-without-spack'], description = 'disable the use of spack environments')
    boolean withoutspack

    @parameter(names=['-offline'], description = 'do not check for remote project updates')
    boolean offline = system.getenv('nxf_offline')=='true'

    @parameter(names=['-entry'], description = 'entry workflow name to be executed', arity = 1)
    string entryname

    @parameter(names=['-main-script'], description = 'the script file to be executed when launching a project directory or repository' )
    string mainscript

    @parameter(names=['-stub-run','-stub'], description = 'execute the workflow replacing process scripts with command stubs')
    boolean stubrun

    @parameter(names=['-preview'], description = "run the workflow script skipping the execution of all processes")
    boolean preview

    @parameter(names=['-latchjit'], description = "generate workflow metadata for a latch execution.")
    boolean latchjit

    @parameter(names=['-plugins'], description = 'specify the plugins to be applied for this run e.g. nf-amazon,nf-tower')
    string plugins

    @parameter(names=['-disable-jobs-cancellation'], description = 'prevent the cancellation of child jobs on execution termination')
    boolean disablejobscancellation

    boolean getdisablejobscancellation() {
        return disablejobscancellation!=null
                ?  disablejobscancellation
                : sysenv.get('nxf_disable_jobs_cancellation') as boolean
    }

    /**
     * optional closure modelling an action to be invoked when the preview mode is enabled
     */
    closure<void> previewaction

    @override
    string getname() { name }

    string getparamsfile() {
        return paramsfile ?: sysenv.get('nxf_params_file')
    }

    boolean hasparams() {
        return params || getparamsfile()
    }

    @override
    void run() {
        final scriptargs = (args?.size()>1 ? args[1..-1] : []) as list<string>
        final pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new abortoperationexception("no project name was specified")

        if( withpodman && withoutpodman )
            throw new abortoperationexception("command line options `-with-podman` and `-without-podman` cannot be specified at the same time")

        if( withdocker && withoutdocker )
            throw new abortoperationexception("command line options `-with-docker` and `-without-docker` cannot be specified at the same time")

        if( withconda && withoutconda )
            throw new abortoperationexception("command line options `-with-conda` and `-without-conda` cannot be specified at the same time")

        if( withspack && withoutspack )
            throw new abortoperationexception("command line options `-with-spack` and `-without-spack` cannot be specified at the same time")

        if( offline && latest )
            throw new abortoperationexception("command line options `-latest` and `-offline` cannot be specified at the same time")

        checkrunname()

        log.info "n e x t f l o w  ~  version ${const.app_ver}"
        plugins.init()

        // -- specify the arguments
        final scriptfile = getscriptfile(pipeline)

        // create the config object
        final builder = new configbuilder()
                .setoptions(launcher.options)
                .setcmdrun(this)
                .setbasedir(scriptfile.parent)
        final config = builder .build()

        // check dsl syntax in the config
        launchinfo(config, scriptfile)

        // check if nxf_ variables are set in nextflow.config
        checkconfigenv(config)

        // -- load plugins
        final cfg = plugins ? [plugins: plugins.tokenize(',')] : config
        plugins.load(cfg)

        // -- load secret provider
        if( secretsloader.isenabled() ) {
            final provider = secretsloader.instance.load()
            config.withsecretprovider(provider)
        }

        // -- create a new runner instance
        final runner = new scriptrunner(config)
        runner.setscript(scriptfile)
        runner.setpreview(this.preview, previewaction)
        runner.setlatchjit(this.latchjit)
        runner.session.profile = profile
        runner.session.commandline = launcher.clistring
        runner.session.ansilog = launcher.options.ansilog
        runner.session.debug = launcher.options.remotedebug
        runner.session.disablejobscancellation = getdisablejobscancellation()

        final istowerenabled = config.navigate('tower.enabled') as boolean
        if( istowerenabled || log.istraceenabled() )
            runner.session.resolvedconfig = configbuilder.resolveconfig(scriptfile.parent, this)
        // note config files are collected during the build process
        // this line should be after `configbuilder#build`
        runner.session.configfiles = builder.parsedconfigfiles
        // set the commit id (if any)
        runner.session.commitid = scriptfile.commitid
        if( this.test ) {
            runner.test(this.test, scriptargs)
            return
        }

        def info = cmdinfo.status( log.istraceenabled() )
        log.debug( '\n'+info )

        // -- add this run to the local history
        runner.verifyandtrackhistory(launcher.clistring, runname)

        // -- run it!
        runner.execute(scriptargs, this.entryname)
    }

    protected checkconfigenv(configmap config) {
        // warn about setting nxf_ environment variables within env config scope
        final env = config.env as map<string, string>
        for( string name : env.keyset() ) {
            if( name.startswith('nxf_') && name!='nxf_debug' ) {
                final msg = "nextflow variables must be defined in the launching environment - the following variable set in the config file is going to be ignored: '$name'"
                log.warn(msg)
            }
        }
    }

    protected void launchinfo(configmap config, scriptfile scriptfile) {
        // -- determine strict mode
        final defstrict = sysenv.get('nxf_enable_strict') ?: false
        final strictmode = config.navigate('nextflow.enable.strict', defstrict)
        if( strictmode ) {
            log.debug "enabling nextflow strict mode"
            nextflowmeta.instance.strictmode(true)
        }
        // -- determine dsl mode
        final dsl = detectdslmode(config, scriptfile.main.text, sysenv)
        nextflowmeta.instance.enabledsl(dsl)
        // -- show launch info
        final ver = nf.dsl2 ? dsl2 : dsl1
        final repo = scriptfile.repository ?: scriptfile.source
        final head = preview ? "* preview * $scriptfile.repository" : "launching `$repo`"
        if( scriptfile.repository )
            log.info "${head} [$runname] dsl${ver} - revision: ${scriptfile.revisioninfo}"
        else
            log.info "${head} [$runname] dsl${ver} - revision: ${scriptfile.getscriptid()?.substring(0,10)}"
    }

    static string detectdslmode(configmap config, string scripttext, map sysenv) {
        // -- try determine dsl version from config file

        final dsl = config.navigate('nextflow.enable.dsl') as string

        // -- script can still override the dsl version
        final scriptdsl = nextflowmeta.checkdslmode(scripttext)
        if( scriptdsl ) {
            log.debug("applied dsl=$scriptdsl from script declararion")
            return scriptdsl
        }
        else if( dsl ) {
            log.debug("applied dsl=$dsl from config declaration")
            return dsl
        }
        // -- if still unknown try probing for dsl1
        if( nextflowmeta.probedsl1(scripttext) ) {
            log.debug "applied dsl=1 by probing script field"
            return dsl1
        }

        final envdsl = sysenv.get('nxf_default_dsl')
        if( envdsl ) {
            log.debug "applied dsl=$envdsl from nxf_default_dsl variable"
            return envdsl
        }
        else {
            log.debug "applied dsl=2 by global default"
            return dsl2
        }
    }

    protected void checkrunname() {
        if( runname == 'last' )
            throw new abortoperationexception("not a valid run name: `last`")
        if( runname && !matchrunname(runname) )
            throw new abortoperationexception("not a valid run name: `$runname` -- it must match the pattern $run_name_pattern")

        if( !runname ) {
            if( historyfile.disabled() )
                throw new abortoperationexception("missing workflow run name")
            // -- make sure the generated name does not exist already
            runname = historyfile.default.generatenextname()
        }

        else if( !historyfile.disabled() && historyfile.default.checkexistsbyname(runname) )
            throw new abortoperationexception("run name `$runname` has been already used -- specify a different one")
    }

    static protected boolean matchrunname(string name) {
        run_name_pattern.matcher(name).matches()
    }

    protected scriptfile getscriptfile(string pipelinename) {
        try {
            getscriptfile0(pipelinename)
        }
        catch (illegalargumentexception | abortoperationexception e) {
            if( e.message.startswith("not a valid project name:") && !guessisrepo(pipelinename)) {
                throw new abortoperationexception("cannot find script file: $pipelinename")
            }
            else
                throw e
        }
    }

    static protected boolean guessisrepo(string name) {
        if( filehelper.geturlprotocol(name) != null )
            return true
        if( name.startswith('/') )
            return false
        if( name.startswith('./') || name.startswith('../') )
            return false
        if( name.endswith('.nf') )
            return false
        if( name.count('/') != 1 )
            return false
        return true
    }

    protected scriptfile getscriptfile0(string pipelinename) {
        assert pipelinename

        /*
         * read from the stdin
         */
        if( pipelinename == '-' ) {
            def file = tryreadfromstdin()
            if( !file )
                throw new abortoperationexception("cannot access `stdin` stream")

            if( revision )
                throw new abortoperationexception("revision option cannot be used when running a script from stdin")

            return new scriptfile(file)
        }

        /*
         * look for a file with the specified pipeline name
         */
        def script = new file(pipelinename)
        if( script.isdirectory()  ) {
            script = mainscript ? new file(mainscript) : new assetmanager().setlocalpath(script).getmainscriptfile()
        }

        if( script.exists() ) {
            if( revision )
                throw new abortoperationexception("revision option cannot be used when running a local script")
            return new scriptfile(script)
        }

        /*
         * try to look for a pipeline in the repository
         */
        def manager = new assetmanager(pipelinename, this)
        def repo = manager.getproject()

        boolean checkforupdate = true
        if( !manager.isrunnable() || latest ) {
            if( offline )
                throw new abortoperationexception("unknown project `$repo` -- note: automatic download from remote repositories is disabled")
            log.info "pulling $repo ..."
            def result = manager.download(revision,deep)
            if( result )
                log.info " $result"
            checkforupdate = false
        }
        // checkout requested revision
        try {
            manager.checkout(revision)
            manager.updatemodules()
            final scriptfile = manager.getscriptfile(mainscript)
            if( checkforupdate && !offline )
                manager.checkremotestatus(scriptfile.revisioninfo)
            // return the script file
            return scriptfile
        }
        catch( abortoperationexception e ) {
            throw e
        }
        catch( exception e ) {
            throw new abortoperationexception("unknown error accessing project `$repo` -- repository may be corrupted: ${manager.localpath}", e)
        }

    }

    static protected file tryreadfromstdin() {
        if( !system.in.available() )
            return null

        getscriptfromstream(system.in)
    }

    static protected file getscriptfromstream( inputstream input, string name = 'nextflow' ) {
        input != null
        file result = file.createtempfile(name, null)
        result.deleteonexit()
        input.withreader { reader reader -> result << reader }
        return result
    }

    @memoized  // <-- avoid parse multiple times the same file and params
    map parsedparams(map configvars) {

        final result = [:]
        final file = getparamsfile()
        if( file ) {
            def path = validateparamsfile(file)
            def type = path.extension.tolowercase() ?: null
            if( type == 'json' )
                readjsonfile(path, configvars, result)
            else if( type == 'yml' || type == 'yaml' )
                readyamlfile(path, configvars, result)
        }

        // set the cli params
        if( !params )
            return result

        for( map.entry<string,string> entry : params ) {
            addparam( result, entry.key, entry.value )
        }
        return result
    }


    static final private pattern dot_escaped = ~/\\\./
    static final private pattern dot_not_escaped = ~/(?<!\\)\./

    static protected void addparam(map params, string key, string value, list path=[], string fullkey=null) {
        if( !fullkey )
            fullkey = key
        final m = dot_not_escaped.matcher(key)
        if( m.find() ) {
            final p = m.start()
            final root = key.substring(0, p)
            if( !root ) throw new abortoperationexception("invalid parameter name: $fullkey")
            path.add(root)
            def nested = params.get(root)
            if( nested == null ) {
                nested = new linkedhashmap<>()
                params.put(root, nested)
            }
            else if( nested !instanceof map ) {
                log.warn "command line parameter --${path.join('.')} is overwritten by --${fullkey}"
                nested = new linkedhashmap<>()
                params.put(root, nested)
            }
            addparam((map)nested, key.substring(p+1), value, path, fullkey)
        }
        else {
            params.put(key.replaceall(dot_escaped,'.'), parseparamvalue(value))
        }
    }


    static protected parseparamvalue(string str) {
        if ( sysenv.get('nxf_disable_params_type_detection') )
            return str

        if ( str == null ) return null

        if ( str.tolowercase() == 'true') return boolean.true
        if ( str.tolowercase() == 'false' ) return boolean.false

        if ( str==~/-?\d+(\.\d+)?/ && str.isinteger() ) return str.tointeger()
        if ( str==~/-?\d+(\.\d+)?/ && str.islong() ) return str.tolong()
        if ( str==~/-?\d+(\.\d+)?/ && str.isdouble() ) return str.todouble()

        return str
    }

    private path validateparamsfile(string file) {

        def result = filehelper.aspath(file)
        def ext = result.getextension()
        if( !valid_params_file.contains(ext) )
            throw new abortoperationexception("not a valid params file extension: $file -- it must be one of the following: ${valid_params_file.join(',')}")

        return result
    }

    static private pattern params_var = ~/(?m)\$\{(\p{javajavaidentifierstart}\p{javajavaidentifierpart}*)}/

    protected string replacevars0(string content, map binding) {
        content.replaceall(params_var) { list<string> matcher ->
            // - the regex matcher is represented as list
            // - the first element is the matching string ie. `${something}`
            // - the second element is the group content ie. `something`
            // - make sure the regex contains at least a group otherwise the closure
            // parameter is a string instead of a list of the call fail
            final placeholder = matcher.get(0)
            final key = matcher.get(1)

            if( !binding.containskey(key) ) {
                final msg = "missing params file variable: $placeholder"
                if(nf.strictmode)
                    throw new abortoperationexception(msg)
                log.warn msg
                return placeholder
            }

            return binding.get(key)
        }
    }

    private void readjsonfile(path file, map configvars, map result) {
        try {
            def text = configvars ? replacevars0(file.text, configvars) : file.text
            def json = (map)new jsonslurper().parsetext(text)
            result.putall(json)
        }
        catch (nosuchfileexception | filenotfoundexception e) {
            throw new abortoperationexception("specified params file does not exists: ${file.touristring()}")
        }
        catch( exception e ) {
            throw new abortoperationexception("cannot parse params file: ${file.touristring()} - cause: ${e.message}", e)
        }
    }

    private void readyamlfile(path file, map configvars, map result) {
        try {
            def text = configvars ? replacevars0(file.text, configvars) : file.text
            def yaml = (map)new yaml().load(text)
            result.putall(yaml)
        }
        catch (nosuchfileexception | filenotfoundexception e) {
            throw new abortoperationexception("specified params file does not exists: ${file.touristring()}")
        }
        catch( exception e ) {
            throw new abortoperationexception("cannot parse params file: ${file.touristring()}", e)
        }
    }

}
