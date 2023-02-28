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

package nextflow.cloud.azure.fusion

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.SysEnv
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionEnv
import org.pf4j.Extension

/**
 * Implements {@link FusionEnv} for AZ cloud
 *
 * @author Pablo Aledo <pablo.aledo@seqera.io>
 */
@Extension
@CompileStatic
class AzFusionEnv implements FusionEnv {

    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {

        final result = new HashMap<String,String>()
        final creds = config.exportAzAccessKeys() ? Global.getAzCredentials() : Collections.<String>emptyList()
        if( creds ) {
            result.AZURE_STORAGE_KEY = creds[0]
        }

        result.AZURE_STORAGE_ACCOUNT = SysEnv.get()['AZURE_STORAGE_ACCOUNT']

        return result
    }
}
