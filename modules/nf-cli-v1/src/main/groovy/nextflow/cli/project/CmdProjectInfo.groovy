/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli.project

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import org.yaml.snakeyaml.Yaml

/**
 * Project info subcommand -- print the information of a downloaded project
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Print project information")
class CmdProjectInfo extends CmdBase {

    static final public String NAME = 'info'

    protected PrintStream out = System.out

    @Parameter(required = true, description = 'project name')
    List<String> args

    @Parameter(names='-d',description = 'Show detailed information', arity = 0)
    boolean detailed

    @Parameter(names='-dd', hidden = true, arity = 0)
    boolean moreDetailed

    @Parameter(names='-o', description = 'Output format, either: text (default), json, yaml')
    String format

    @Parameter(names=['-u','-check-updates'], description = 'Check for remote updates')
    boolean checkForUpdates

    @Override
    String getName() { NAME }

    @Override
    void run() {
        final level = moreDetailed ? 2 : ( detailed ? 1 : 0 )
        Plugins.init()
        try (final manager = new AssetManager(args[0])) {
            if( manager.isNotInitialized() ) {
                throw new AbortOperationException("Unknown project `${args[0]}`")
            }

            if( !format || format == 'text' ) {
                printText(manager,level)
                return
            }

            final map = createMap(manager)
            if( format == 'json' ) {
                printJson(map)
            }
            else if( format == 'yaml' ) {
                printYaml(map)
            }
            else
                throw new AbortOperationException("Unknown output format: $format");
        }
    }

    protected printText(AssetManager manager, int level) {
        final manifest = manager.getManifest()

        out.println " project name: ${manager.project}"
        out.println " repository  : ${manager.repositoryUrl}"
        out.println " local path  : ${manager.projectPath}"
        out.println " main script : ${manager.mainScriptName}"
        if( manager.homePage && manager.homePage != manager.repositoryUrl )
            out.println " home page   : ${manager.homePage}"
        if( manifest.description )
            out.println " description : ${manifest.description}"
        if( manifest.author )
            out.println " author      : ${manifest.author}"

        def revs = manager.getRevisions(level)
        if( revs.size() == 1 )
            out.println " revision    : ${revs[0]}"
        else {
            out.println " revisions   : "
            revs.each { out.println " $it" }
        }

        if( !checkForUpdates )
            return

        def updates = manager.getUpdates(level)
        if( updates ) {
            if( updates.size() == 1 && revs.size() == 1 )
                out.println " updates     : ${updates[0]}"
            else {
                out.println " updates     : "
                updates.each { out.println " $it" }
            }
        }

        out.flush()
    }

    protected Map createMap(AssetManager manager) {
        def result = [:]
        result.projectName = manager.project
        result.repository = manager.repositoryUrl
        result.localPath = manager.projectPath?.toString()
        result.manifest = manager.manifest.toMap()
        result.revisions = manager.getBranchesAndTags(checkForUpdates)
        return result
    }

    protected printJson(Map map) {
        out.println JsonOutput.prettyPrint(JsonOutput.toJson(map))
        out.flush()
    }

    protected printYaml(Map map) {
        out.println new Yaml().dump(map).toString()
        out.flush()
    }
}
