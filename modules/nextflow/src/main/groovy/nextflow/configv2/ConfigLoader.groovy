/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.configv2

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.micronaut.context.ApplicationContext
import io.micronaut.context.env.Environment

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigLoader {

    private ApplicationContext context

    ConfigLoader load(List<Path> configs) {
        final configPaths = configs.collect(it-> 'file:'+it.toAbsolutePath()).join(',')
        log.debug "Loading config files: ${configPaths}"
        System.setProperty('micronaut.config.files', "${configPaths}")

        final props = ['micronaut.bootstrap.context': 'false'] as Map<String,Object>
        context = ApplicationContext.run(props, Environment.CLI)

        return this
    }

    NextflowOpts nextflowConfig() {
        context.getBean(NextflowOpts)
    }

}
