/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cli

import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import org.fusesource.jansi.Ansi

/**
 * Interface for top-level CLI options.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
trait ILauncherOptions {

    abstract Boolean getAnsiLogCli()

    abstract boolean isBackground()

    abstract List<String> getConfig()

    abstract List<String> getDebug()

    abstract boolean getIgnoreConfigIncludes()

    abstract String getLogFile()

    abstract boolean isQuiet()

    abstract String getSyslog()

    abstract List<String> getTrace()

    abstract List<String> getUserConfig()

    abstract void setAnsiLog(boolean value)

    abstract void setBackground(boolean value)

    boolean getAnsiLog() {
        if( ansiLogCli && quiet )
            throw new AbortOperationException("Command line options `quiet` and `ansi-log` cannot be used together")

        if( ansiLogCli != null )
            return ansiLogCli

        if( background )
            return ansiLogCli = false

        if( quiet )
            return ansiLogCli = false

        final env = System.getenv('NXF_ANSI_LOG')
        if( env ) try {
            return Boolean.parseBoolean(env)
        }
        catch (Exception e) {
            log.warn "Invalid boolean value for variable NXF_ANSI_LOG: $env -- it must be 'true' or 'false'"
        }
        return Ansi.isEnabled()
    }

    boolean hasAnsiLogFlag() {
        ansiLogCli==true || System.getenv('NXF_ANSI_LOG')=='true'
    }

}
