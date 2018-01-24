/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import nextflow.extension.Bolts

/**
 * Validation helper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CheckHelper {

    /**
     * Valid a method named parameters map
     *
     * @param name The name of the method, only in the error reported message
     * @param params The actual parameters map
     * @param valid A map providing for each parameter name the valid values
     * @throws IllegalArgumentException when the parameter include an unexpected parameter name or value
     */
    static void checkParams( String name, Map<String,?> params, Map<String,?> valid )  {

        if( !params ) return

        def allKeys = valid.keySet()
        for( String key : params.keySet() ) {
            if( !allKeys.contains(key) )
                throw new IllegalArgumentException("Unknown argument '${key}' for operator '$name' -- Possible arguments: ${allKeys.join(', ')}")

            final value = params.get(key)
            final accepted = valid.get(key)
            if( accepted instanceof Collection ) {
                boolean ok = false
                final itr = accepted.iterator()
                while( !ok && itr.hasNext() )
                    ok |= isValid(value, itr.next())
                if( !ok )
                    throw new IllegalArgumentException("Value '${value}' cannot be used in in parameter '${key}' for operator '$name' -- Possible values: ${(accepted as Collection).join(', ')}")
            }

            else if( !isValid(value, accepted) )
                throw new IllegalArgumentException("Value '${value}' cannot be used in in parameter '${key}' for operator '$name' -- Value don't match: ${accepted}")

        }

    }

    /**
     * Check if the provide value is included in the specified range
     *
     * @param value A value to verify
     * @param range The range it may be a {@link Class} a {@link Collection} of values, a regexp {@link Pattern} or a specific value
     * @return {@code true} is the match is satisfied or {@code false} otherwise
     */
    static boolean isValid( value, range ) {
        if( range instanceof Class )
            return range.isAssignableFrom(value.class)

        if( range instanceof Collection )
            return range.contains(value)

        if( value != null && range instanceof Pattern )
            return range.matcher(value.toString()).matches()

        value == range
    }

    /**
     *  Verify that all method named parameters are included in the provided list
     *
     * @param name The method name used in the reported error message
     * @param params The list of accepted named parameters
     * @param valid The list of accepted parameter names
     * @throws IllegalArgumentException a parameter is not included the in valid names list
     */
    static void checkParams( String name, Map<String,?> params, List<String> valid )  {
        if( !params ) return

        for( String key : params.keySet() ) {
            if( !valid.contains(key) ) {
                def matches = Bolts.closest(valid,key) ?: valid
                def message = "Unknown argument '${key}' for operator '$name'. Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
                throw new IllegalArgumentException(message)
            }
        }
    }


    /**
     *  Verify that all method named parameters are included in the provided list
     *
     * @param name The method name used in the reported error message
     * @param params The list of accepted named parameters
     * @param valid The list of accepted parameter names
     * @throws IllegalArgumentException a parameter is not included the in valid names list
     */
    static void checkParams( String name, Map<String,?> params, String... valid ) {
        checkParams(name, params, Arrays.asList(valid))
    }


}
