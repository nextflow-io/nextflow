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
        final list = Plugins.getExtensions(FusionEnv)
        final result = new HashMap<String,String>()
        final env = getPluginsEnv(list, scheme, config)
        if( env )
            result.putAll(env)
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

    protected Map<String,String> getPluginsEnv(List<FusionEnv> list, String scheme, FusionConfig config) {
        try {
            return getPluginsEnv0(list,scheme,config)
        }
        catch (ReportWarningException e) {
            log.warn1(e.message, causedBy:e, cacheKey:e.kind)
            return Map.of()
        }
    }

    protected Map<String,String> getPluginsEnv0(List<FusionEnv> list, String scheme, FusionConfig config) {
        final result = new HashMap()
        for( FusionEnv it : list ) {
            final env = it.getEnvironment(scheme,config)
            if( env )
                result.putAll(env)
        }
        if( !result.containsKey('FUSION_LICENSE_TOKEN') ) {
            final msg = "Unable to validate Fusion license - Make sure to have added in your config 'tower.accessToken' or have defined the variable TOWER_ACCESS_TOKEN in your launch environment"
            throw new ReportWarningException(msg, 'getFusionLicenseException')
        }
        return result
    }

}
