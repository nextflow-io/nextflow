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

package nextflow.config

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.scm.ProviderConfig
import nextflow.util.StringUtils

/**
 * Reads and parses the SCM config file, which needs the runtime config parser
 * and therefore cannot live with {@link ProviderConfig} in the nf-scm module
 */
@Slf4j
@CompileStatic
class ProviderConfigParser {

    static Map getDefault() {
        return getFromFile(ProviderConfig.getScmConfigPath())
    }

    static Map getFromFile(Path file) {
        try {
            // note: since this can be a remote file via e.g. read via HTTP
            // it should not read more than once
            final content = file.text
            dumpScmContent(file, content)
            final result = parse(content)
            dumpConfig(result)
            return result
        }
        catch (NoSuchFileException | FileNotFoundException e) {
            if( file == ProviderConfig.DEFAULT_SCM_FILE ) {
                return new LinkedHashMap()
            }
            else
                throw new AbortOperationException("Missing SCM config file: ${file.toUriString()} - Check the env variable NXF_SCM_FILE")
        }
        catch ( UnknownHostException e ) {
            final message = "Unable to access config file '${file?.toUriString()}' -- Unknown host: ${e}"
            throw new ConfigParseException(message,e)
        }
        catch( IOException e ) {
            final message = "Unable to access config file '${file?.toUriString()}' -- Cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
        catch( Exception e ) {
            final message = "Failed to parse config file '${file?.toUriString()}' -- Cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
    }

    static Map parse(String text) {
        def slurper = ConfigParserFactory.create()
        slurper.setBinding(ProviderConfig.env)
        return slurper.parse(text)
    }

    static private void dumpScmContent(Path file, String content) {
        try {
            log.trace "Parsing SCM config path: ${file.toUriString()}\n${StringUtils.stripSecrets(content)}\n"
        }catch(Exception e){
            log.debug "Error dumping configuration ${file.toUriString()}", e
        }
    }

    static private void dumpConfig(Map config) {
        try {
            log.debug "Detected SCM config: ${StringUtils.stripSecrets(config).toMapString()}"
        }
        catch (Throwable e) {
            log.debug "Failed to dump SCM config: ${e.message ?: e}", e
        }
    }
}
