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

package nextflow.sql.config

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model SQL dataSources configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode(includeFields = true)
@CompileStatic
class SqlConfig {

    private Map<String,SqlDataSource> dataSources
    private SqlDataSource defSource

    SqlConfig(Map config) {
        dataSources = parseDataSources(config)
    }

    SqlDataSource getDataSource(String name) {
        return name=='default'
                ? defSource
                : dataSources.get(name)
    }

    List<String> getDataSourceNames() {
        return new ArrayList<String>(dataSources.keySet())
    }

    protected Map parseDataSources(Map<String,Map> config) {
        // create the `default`datasource as fallback for missing values
        this.defSource = new SqlDataSource(config?.'default' as Map ?: Collections.emptyMap())

        // setup other db by name
        final result = new LinkedHashMap<String,SqlDataSource>()
        for( Map.Entry<String,Map> entry : config ) {
            if( entry.key == 'default' ) {
                result.put( entry.key, defSource )
            }
            else {
                result.put( entry.key, new SqlDataSource(entry.value, defSource) )
            }
        }

        return result
    }

}
