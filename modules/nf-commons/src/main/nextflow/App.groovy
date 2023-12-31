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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.micronaut.context.ApplicationContext
import io.micronaut.context.env.Environment
import nextflow.plugin.PluginService
/**
 * Creates Micronaut application context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class App {

    private static Object lock = new Object()

    private static App app

    private ApplicationContext context

    static App getInstance() {
        if( !app )
            throw new IllegalStateException("Application context has not been started")
        return app
    }

    static App start() {
        synchronized (lock) {
            if( app )
                throw new IllegalStateException("Application context has been already started")
            return app = new App()
        }
    }

    static void shutdown(boolean quietly=false) {
        if( !app ) {
            final msg = "Application context has not been started"
            if( quietly ) {
                log.debug msg
                return
            }
            throw new IllegalStateException(msg)
        }
        app.context.stop()
        app = null
    }

    private App() {
        final props = ['micronaut.bootstrap.context': 'false'] as Map<String,Object>
        context = ApplicationContext.run(props, Environment.CLI)
    }

    ApplicationContext context() {
        return context
    }

    <T> T getBean(Class<T> clazz) {
        context.getBean(clazz)
    }

    PluginService getPluginService() {
        context.getBean(PluginService)
    }

}
