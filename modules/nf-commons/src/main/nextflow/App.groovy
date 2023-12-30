/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow

import groovy.transform.Memoized
import io.micronaut.context.ApplicationContext
import io.micronaut.context.env.Environment
import nextflow.plugin.PluginService

/**
 * Creates Micronaut application context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class App {

    @Memoized
    static synchronized ApplicationContext context() {
        final props = ['micronaut.bootstrap.context': 'false'] as Map<String,Object>
        final result = ApplicationContext.run(props, Environment.CLI)
        return result
    }

    static <T> T get(Class<T> clazz) {
        context().getBean(clazz)
    }

    static PluginService getPluginService() {
        get(PluginService)
    }
}
