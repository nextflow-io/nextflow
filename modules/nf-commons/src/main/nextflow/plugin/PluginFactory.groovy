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

package nextflow.plugin

import io.micronaut.context.annotation.Context
import io.micronaut.context.annotation.Factory
import io.micronaut.context.annotation.Primary
/**
 * Implements a factory class for {@link PluginService}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Factory
class PluginFactory {

    @Primary
    @Context
    PluginService createStandardService() {
        final result = new PluginsFacade()
        result.init(false)
        return result
    }


//    @Named('embeddable')
//    @Singleton
//    PluginService createEmbeddableService() {
//        final result = new PluginsFacade()
//        result.init(true)
//        return result
//    }
}
