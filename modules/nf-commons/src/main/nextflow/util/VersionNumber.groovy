/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Model a semantic version number
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class VersionNumber implements Comparable {

    static private Pattern CHECK = ~/\s*([!<=>]+)?\s*([0-9A-Za-z\-_\.]+)(\+)?/

    private List<String> version

    /**
     * Create a version number object
     *
     * @param str A dot separated version number string eg {@code 2.3.12}
     *
     */
    VersionNumber(String str) {
        version = str ? str.tokenize('.') : ['0']
    }

    /**
     * @return The number of components in the version number eg {@code new VersionNumber('2.3.12').size() == 3}
     */
    int size() { version.size() }

    /**
     * Retrieve the i-th component in a version number
     *
     * @param i The i-th index of the component to retrieve. The first is 0.
     * @return The i-th version number component string
     */
    String getAt( int i ) { version[i] }

    /**
     * @return The major version number ie. the first component
     */
    String getMajor() { version[0] }

    /**
     * @return The minor version number ie. the second component
     */
    String getMinor() { version[1] }

    /**
     * @return The minor version number ie. the third component
     */
    String getPatch() { version[2] }

    /**
     * Implements version number comparator logic
     *
     * @param   o the object to be compared.
     * @return  a negative integer, zero, or a positive integer as this object
     *          is less than, equal to, or greater than the specified object.
     */
    @Override
    int compareTo(Object obj) {
        if( obj == null )
            return 1

        if( obj instanceof VersionNumber ) {
            return compare(obj)
        }
        if( obj instanceof CharSequence ) {
            return compare(new VersionNumber(obj.toString()))
        }
        if( obj instanceof Number ) {
            return compare(new VersionNumber(obj.toString()))
        }

        throw new ClassCastException("Cannot compare VersionNumber with object $obj [${obj.class}]")
    }

    private int compare( VersionNumber obj ) {
        def _this = version
        def _that = obj.version
        while ( _this.size() < _that.size() )
            _this << '0'
        while( _this.size() > _that.size() )
            _that << '0'

        for( int i=0; i<_this.size(); i++ ) {
            def result = compare(_this[i], _that[i])
            if( result != 0 )
                return result
        }

        return 0
    }

    private int compare( String a, String b ) {
        return (
                a.isInteger() && b.isInteger()
                    ? a.toInteger() <=> b.toInteger()
                    : a <=> b
                )
    }

    String toString() {
        version.join('.')
    }

    @Override
    int hashCode() {
        Objects.hash( version as Object[] )
    }

    @Override
    boolean equals(Object other) {
        compareTo(other) == 0
    }

    private boolean checkPlus(String checkStr) {
        if( checkStr.endsWith('.') )
            checkStr += '0'
        def check = new VersionNumber(checkStr)
        def last = check.size()-1
        def List<String> _this = this.version
        while ( _this.size() < check.size() )
            _this << '0'

        for( int i=0; i<last; i++ ) {
            if( _this[i] != check[i] ) {
                return false
            }

        }
        def result = compare(check[last], _this[last])
        return result <= 0
    }

    /**
     * Checks if the version number satisfies the specified condition.
     * The condition can be specified by prefix the request version with
     * traditional logical operators eg:
     *
     * <pre>
     *    assert new VersionNumber('1.2').check('< 1.4')
     *    assert new VersionNumber('2.3').check('>= 2.0')
     * </pre>
     *
     * Or by post-fixing the request version by using a `+` sign, eg:
     *
     * <pre>
     *    assert new VersionNumber('1.2').check('1.1+')
     *    assert new VersionNumber('2.3').check('2.+')
     * </pre>
     *
     * @param condition
     * @return
     */
    boolean matches(String condition) {
        for( String token : condition.tokenize(',')) {
            if( !matches0(token) )
                return false
        }
        return true
    }

    private boolean matches0(String condition) {
        def matcher = CHECK.matcher(condition)
        if( !matcher.matches() )
            throw new IllegalArgumentException("Not a valid version check condition: $condition")

        def operator = matcher.group(1)
        def version = matcher.group(2)
        def plus = matcher.group(3)
        if( operator && plus )
            throw new IllegalArgumentException("Not a valid version check condition: $condition")

        if( plus )
            return checkPlus(version)
        if( !operator )
            operator = '='

        switch (operator) {
            case '=':
            case '==':
                return compareTo(version) == 0

            case '<':
                return compareTo(version) < 0

            case '<=':
                return compareTo(version) <= 0

            case '>':
                return compareTo(version) > 0

            case '>=':
                return compareTo(version) >= 0

            case '!=':
            case '<>':
                return compareTo(version) != 0

            default:
                throw new IllegalArgumentException("Not a version check operator: $operator")
        }
    }

}
