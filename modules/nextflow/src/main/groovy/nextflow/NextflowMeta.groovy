/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow

import static nextflow.extension.Bolts.*

import java.text.SimpleDateFormat

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.util.VersionNumber
/**
 * Models nextflow script properties and metadata
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Singleton(strict = false)
@ToString(includeNames = true)
@EqualsAndHashCode
class NextflowMeta {

    static trait Flags {
        abstract float dsl
        abstract boolean strict
        abstract boolean moduleBinaries
    }

    @Slf4j
    static class Preview implements Flags {
        @Deprecated volatile float dsl
        @Deprecated boolean strict
        boolean recursion
        boolean moduleBinaries
        boolean types

        void setRecursion(Boolean recursion) {
            if( recursion )
                log.warn "NEXTFLOW RECURSION IS A PREVIEW FEATURE - SYNTAX AND FUNCTIONALITY CAN CHANGE IN FUTURE RELEASES"
            this.recursion = recursion
        }
    }

    static class Features implements Flags {
        volatile float dsl
        boolean strict
        boolean moduleBinaries
    }

    final VersionNumber version
    final int build

    /*
     * Timestamp as dd-MM-yyyy HH:mm UTC formatted string
     */
    final String timestamp

    final Preview preview = new Preview()

    final Features enable = new Features()

    private NextflowMeta() {
        version = new VersionNumber(BuildInfo.version)
        build = BuildInfo.buildNum as int
        timestamp = BuildInfo.timestampUTC
    }

    protected NextflowMeta(String ver, int build, String timestamp ) {
        this.version = new VersionNumber(ver)
        this.build = build
        this.timestamp = timestamp
    }

    Map featuresMap() {
        final result = new LinkedHashMap()
        if( isStrictModeEnabled() )
            result.strict = true
        if( isModuleBinariesEnabled() )
            result.moduleBinaries = true
        return result
    }

    Map toJsonMap() {
        final result = new LinkedHashMap<>(5)
        result.version = version.toString()
        result.build = build
        result.timestamp = parseDateStr(timestamp)
        result.enable = featuresMap()
        return result
    }

    private Date parseDateStr(String str) {
        def fmt = new SimpleDateFormat(DATETIME_FORMAT + ' Z')
        fmt.parse(str)
    }

    boolean isStrictModeEnabled() {
        return enable.strict
    }

    void strictMode(boolean mode) {
        enable.strict = mode
    }

    boolean isModuleBinariesEnabled() {
        return enable.moduleBinaries
    }

    void moduleBinaries(boolean mode) {
        enable.moduleBinaries = mode
    }

}
