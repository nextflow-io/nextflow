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

package nextflow.fusion

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ReportWarningException
import nextflow.plugin.Plugins
/**
 * Provider strategy for {@link FusionEnv}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FusionEnvProvider {

    Map<String,String> getEnvironment(String scheme) {
        final config = FusionConfig.getConfig()
        final result = new HashMap<String,String>()
        final env = getFusionEnvironment(scheme,config)
        result.putAll(env)
        final validator = getFusionToken(scheme,config)
        result.putAll(validator)
        // tags setting
        if( config.tagsEnabled() )
            result.FUSION_TAGS = config.tagsPattern()
        // logs setting
        if( config.logOutput() )
            result.FUSION_LOG_OUTPUT = config.logOutput()
        if( config.logLevel() )
            result.FUSION_LOG_LEVEL = config.logLevel()
        if( config.cacheSize() )
            result.FUSION_CACHE_SIZE = "${config.cacheSize().toMega()}M"
        return result
    }

    protected Map<String,String> getFusionEnvironment(String scheme, FusionConfig config) {
        final list = Plugins.getExtensions(FusionEnv)
        log.debug "Fusion environment extensions=$list"
        final result = new HashMap<String,String>()
        for( FusionEnv it : list ) {
            final env = it.getEnvironment(scheme,config)
            if( env )
                result.putAll(env)
        }
        return result
    }

    protected Map<String,String> getFusionToken(String scheme, FusionConfig config) {
        final providers = Plugins.getPriorityExtensions(FusionToken)
        log.debug "Fusion token extensions=$providers"
        final result = new HashMap<String,String>()
        try {
            final env = providers.first().getEnvironment(scheme,config)
            if( env )
                result.putAll(env)
        }
        catch (ReportWarningException e) {
            log.warn1(e.message, cacheKey:e.kind)
        }
        return result
    }

}
