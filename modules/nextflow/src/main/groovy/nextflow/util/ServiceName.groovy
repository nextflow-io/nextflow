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

package nextflow.util
import java.lang.annotation.ElementType
import java.lang.annotation.Retention
import java.lang.annotation.RetentionPolicy
import java.lang.annotation.Target
/**
 * Assign a name to a runtime loaded service object.
 *
 * Beside the generic name "service", this class is only mean for
 * annotating Nextflow executor classes.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
@interface ServiceName {

    /**
     * @return The name of the executor service.
     */
    String value()

    /**
     * Determine if the service should have precedence over a service
     * having #important flag marked as false.
     *
     * NOTE: this option is deprecated and it's only meant to be used for
     * Nextflow executor discovery mechanism.
     *
     * Consider using {@link nextflow.plugin.PluginsFacade#getPriorityExtensions}
     * for new services or plugin extensions.
     *
     * @return {@code true} when the executor should have priority over non-important ones
     */
    @Deprecated
    boolean important() default Boolean.FALSE

}
