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
package nextflow.data.cid

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import nextflow.data.config.DataConfig
import nextflow.plugin.Priority

/**
 * Default Factory for CidStore.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Priority(0)
class DefaultCidStoreFactory extends CidStoreFactory {

    private static Pattern SCHEME = ~/^([a-zA-Z][a-zA-Z\d+\-.]*):/
    private static List SUPPORTED_SCHEMES = ['file', 's3', 'gs', 'az']

    @Override
    boolean canOpen(DataConfig config) {
        final loc = config.store.location
        if( !loc ) {
            return true
        }
        final matcher = SCHEME.matcher(loc)
        return matcher.find() ? matcher.group(1) in SUPPORTED_SCHEMES : true
    }

    @Override
    protected CidStore newInstance(DataConfig config) {
        return new DefaultCidStore() .open(config)
    }

}
