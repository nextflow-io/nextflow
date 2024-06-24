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

package nextflow.cli

import static nextflow.scm.AssetManager.DEFAULT_REVISION_DIRNAME

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.scm.AssetManager

/**
 * CLI sub-command LIST. Prints a list of locally installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "List all downloaded projects")
class CmdList extends CmdBase {

    static final public NAME = 'list'

    @Parameter(names=['-a','-all-revisions'], description = 'For each project, also list revisions')
    Boolean allRevisions

    @Parameter(names=['-all-commits'], description = 'For each project, also list all downloaded commits')
    Boolean allCommits

    @Parameter(names='-d',description = 'Show commit information for revisions', arity = 0)
    boolean detailed

    @Parameter(names='-dd', hidden = true, arity = 0)
    boolean moreDetailed

    @Override
    final String getName() { NAME }

    @Override
    void run() {

        def all = AssetManager.list()
        if( !all ) {
            log.info '(none)'
            return
        }

    if( moreDetailed )
        detailed = true
    if( detailed && allRevisions ) {
        all.each{
            println(" $it")
            def revManager = new AssetManager(it)
            revManager.listRevisionsAndCommits().each{ k,v ->
                if( k == DEFAULT_REVISION_DIRNAME )
                    k = '(default)'
                if( !moreDetailed )
                    v = v.substring(0,10)
                println("   $v $k") }
        }
    }
    else if( allRevisions ) {
        all.each{
            println(" $it")
            def revManager = new AssetManager(it)
            revManager.listRevisions().each{
                if( it == DEFAULT_REVISION_DIRNAME )
                    it = '(default)'
                println("   $it")
            }
        }
    } else if( allCommits ) {
        all.each{
            println(" $it")
            def revManager = new AssetManager(it)
            revManager.listCommits().each{ println("   $it") }
        }
    } else {
        all.each{ println(" $it") }
    }
    }

}
