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

import java.nio.file.Paths

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import static PluginExecAware.CMD_SEP

/**
 * Abstract implementation for plugin commands
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait PluginAbstractExec implements PluginExecAware {

    private static Logger log = LoggerFactory.getLogger(PluginAbstractExec)

    private Session session
    private Launcher launcher

    Session getSession() { session }

    Launcher getLauncher() { launcher }

    abstract List<String> getCommands()

    @Override
    final int exec(Launcher launcher1, String pluginId, String cmd, List<String> args) {
        this.launcher = launcher1
        // create the config
        final config = new ConfigBuilder()
                .setOptions(launcher1.options)
                .setBaseDir(Paths.get('.'))
                .build()

        if( !cmd || cmd !in getCommands() ) {
            def msg = cmd ? "Invalid command '$pluginId:$cmd' - usage: nextflow plugin ${pluginId}${CMD_SEP}<command>" : "Usage: nextflow plugin ${pluginId}${CMD_SEP}<command>"
            msg += "\nAvailable commands:"
            for( String it : getCommands() )
                msg += "\n $it"
            throw new AbortOperationException(msg)
        }

        // create the session object
        this.session = new Session(config)
        // invoke the command
        try {
            exec(cmd, args)
        }
        catch (Throwable e) {
            log.error("Unexpected error on command: $cmd - cause: ${e.message}", e)
        }
        finally {
            session.destroy()
        }
        return 0
    }

    abstract int exec(String cmd, List<String> args)


}
