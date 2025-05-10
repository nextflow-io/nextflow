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
 *
 */

package nextflow.cli

import static nextflow.cli.PluginExecAware.CMD_SEP

import java.nio.file.Path

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.plugin.util.PluginRefactor
import org.eclipse.jgit.api.Git
/**
 * Plugin manager command
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Execute plugin-specific commands")
class CmdPlugin extends CmdBase {

    @Override
    String getName() {
        return 'plugin'
    }

    @DynamicParameter(names = "--", description = "Custom plugin parameters go here", hidden = true)
    private Map<String, String> params = new HashMap<>();

    @Parameter(hidden = true)
    List<String> args

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("Missing plugin command - usage: nextflow plugin install <pluginId,..>")
        // setup plugins system
        Plugins.init()
        Runtime.addShutdownHook((it)-> Plugins.stop())
        
        // check for the plugins install
        if( args[0] == 'install' ) {
            if( args.size()!=2 )
                throw new AbortOperationException("Missing plugin install target - usage: nextflow plugin install <pluginId,..>")
            Plugins.pull(args[1].tokenize(','))
        }
        else if( args[0] == 'create' ) {
            createPlugin(args)
        }
        // plugin run command
        else if( args[0].contains(CMD_SEP) ) {
            final head = args.pop()
            final items = head.tokenize(CMD_SEP)
            final target = items[0]
            final cmd = items[1] ? items[1..-1].join(CMD_SEP) : null

            // push back the command as the first item
            Plugins.start(target)
            final wrapper = Plugins.manager.getPlugin(target)
            if( !wrapper )
                throw new AbortOperationException("Cannot find target plugin: $target")
            final plugin = wrapper.getPlugin()
            if( plugin instanceof PluginExecAware ) {
                def mapped = [] as List<String>
                params.entrySet().each{
                    mapped << "--$it.key".toString()
                    mapped << "$it.value".toString()
                }
                args.addAll(mapped)
                final ret = plugin.exec(getLauncher(), target, cmd, args)
                // use explicit exit to invoke the system shutdown hooks
                System.exit(ret)
            }
            else
                throw new AbortOperationException("Invalid target plugin: $target")
        }
        else {
            throw new AbortOperationException("Invalid plugin command: ${args[0]}")
        }
    }

    static createPlugin(List<String> args) {
        if( args != ['create'] && (args[0] != 'create' || !(args.size() in [3, 4])) )
            throw new AbortOperationException("Invalid create parameters - usage: nextflow plugin create <Plugin name> <Organization name>")

        final refactor = new PluginRefactor()
        if( args.size()>1 ) {
            refactor.withPluginName(args[1])
            refactor.withOrgName(args[2])
            refactor.withPluginDir(Path.of(args[3] ?: refactor.pluginName).toFile())
        }
        else {
            // Prompt for plugin name
            print "Enter plugin name: "
            refactor.withPluginName(readLine())

            // Prompt for maintainer organization
            print "Enter organization: "

            // Prompt for plugin path (default to the normalised plugin name)
            refactor.withOrgName(readLine())
            print "Enter project path [${refactor.pluginName}]: "
            refactor.withPluginDir(Path.of(readLine() ?: refactor.pluginName).toFile())

            // confirm and proceed
            print "All good, are you OK to continue [y/N]? "
            final confirm = readLine()
            if( confirm!='y' )
                return
        }

        // the final directory where the plugin is created
        final File targetDir = refactor.getPluginDir()

        // clone the template repo
        clonePluginTemplate(targetDir)
        // now refactor the template code
        refactor.apply()
        // remove git plat
        cleanup(targetDir)
        // done
        println "Plugin created successfully at path: $targetDir"
    }

    static private String readLine() {
        final console = System.console()
        return console != null
            ? console.readLine()
            : new BufferedReader(new InputStreamReader(System.in)).readLine()
    }

    static private void clonePluginTemplate(File targetDir) {
        final templateUri = "https://github.com/nextflow-io/nf-plugin-template.git"
        try {
            Git.cloneRepository()
                .setURI(templateUri)
                .setDirectory(targetDir)
                .setBranchesToClone(["refs/tags/v0.1.0"])
                .setBranch("refs/tags/v0.1.0")
                .call()
        }
        catch (Exception e) {
            throw new AbortOperationException("Unable to clone pluging template repository - cause: ${e.message}")
        }
    }

    static private void cleanup(File targetDir) {
        new File(targetDir, '.git').deleteDir()
        new File(targetDir, '.github').deleteDir()
    }
}
