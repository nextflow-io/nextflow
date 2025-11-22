/*
 * Copyright 2023, Seqera Labs
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

package nextflow.script.params

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.exception.IllegalArityException

/**
 * Implements an arity option for process inputs and outputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
trait ArityParam {

    Range arity

    Range getArity() { arity }

    def setArity(String value) {
        if( value.isInteger() ) {
            def n = value.toInteger()
            this.arity = new Range(n, n)
            return this
        }

        final tokens = value.tokenize('..')
        if( tokens.size() == 2 ) {
            final min = tokens[0]
            final max = tokens[1]
            if( min.isInteger() && (max == '*' || max.isInteger()) ) {
                this.arity = new Range(
                    min.toInteger(),
                    max == '*' ? Integer.MAX_VALUE : max.toInteger()
                )
                return this
            }
        }

        throw new IllegalArityException("Path arity should be a number (e.g. '1') or a range (e.g. '1..*')")
    }

    /**
     * Determine whether a single output file should be unwrapped.
     */
    boolean isSingle() {
        return !arity || arity.max == 1
    }

    boolean isValidArity(int size) {
        return !arity || arity.contains(size)
    }

    @EqualsAndHashCode
    static class Range {
        int min
        int max

        Range(int min, int max) {
            if( min<0 )
                throw new IllegalArityException("Path arity min value must be greater or equals to 0")
            if( max<1 )
                throw new IllegalArityException("Path arity max value must be greater or equals to 1")
            if( min==0 && max==1 )
                throw new IllegalArityException("Path arity 0..1 is not allowed")
            this.min = min
            this.max = max
        }

        boolean contains(int value) {
            min <= value && value <= max
        }

        @Override
        String toString() {
            min == max
                ? min.toString()
                : "${min}..${max == Integer.MAX_VALUE ? '*' : max}".toString()
        }
    }

}
