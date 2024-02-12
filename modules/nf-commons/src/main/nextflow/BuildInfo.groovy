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

package nextflow

import static nextflow.extension.Bolts.*

import java.text.SimpleDateFormat

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Model app build information
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BuildInfo {

    private static Properties properties

    static {
        final BUILD_INFO = '/META-INF/build-info.properties'
        properties = new Properties()
        try {
            properties.load( BuildInfo.getResourceAsStream(BUILD_INFO) )
        }
        catch( Exception e ) {
            log.warn "Unable to parse $BUILD_INFO - Cause ${e.message ?: e}"
        }
    }

    static Properties getProperties() { properties }

    static String getVersion() { properties.getProperty('version') }

    static String getCommitId() { properties.getProperty('commitId')}

    static String getBuildNum() { properties.getProperty('build') }

    static long getTimestampMillis() { properties.getProperty('timestamp') as long }

    static String getTimestampUTC() {
        def tz = TimeZone.getTimeZone('UTC')
        def fmt = new SimpleDateFormat(DATETIME_FORMAT)
        fmt.setTimeZone(tz)
        fmt.format(new Date(getTimestampMillis())) + ' ' + tz.getDisplayName( true, TimeZone.SHORT )
    }

    static String getTimestampLocal() {
        def tz = TimeZone.getDefault()
        def fmt = new SimpleDateFormat(DATETIME_FORMAT)
        fmt.setTimeZone(tz)
        fmt.format(new Date(getTimestampMillis())) + ' ' + tz.getDisplayName( true, TimeZone.SHORT )
    }

    static String getTimestampDelta() {
        if( getTimestampUTC() == getTimestampLocal() ) {
            return ''
        }

        final utc = getTimestampUTC().tokenize(' ')
        final loc = getTimestampLocal().tokenize(' ')

        final result = utc[0] == loc[0] ? loc[1,-1].join(' ') : loc.join(' ')
        return "($result)"
    }

    static String getFullVersion() {
        "${version}_${commitId}"
    }

}
