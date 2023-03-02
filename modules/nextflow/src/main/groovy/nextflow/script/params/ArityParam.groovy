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


/**
 * Implements an arity option for process inputs and outputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
trait ArityParam {

    private Range arity = null

    def setArity(String arity) {
        if( arity.isInteger() ) {
            def n = arity.toInteger()
            this.arity = new Range(n, n)
        }
        else if( arity.contains('..') && arity.tokenize('..').size() == 2 ) {
            def (min, max) = arity.tokenize('..')
            this.arity = new Range(
                min == '*' ? 0 : min.toInteger(),
                max == '*' ? Integer.MAX_VALUE : max.toInteger()
            )
        }
        else {
            throw new IllegalArgumentException("Output path arity should be a number (e.g. '1') or range (e.g. '1..*')")
        }
        return this
    }

    Range getArity() {
        arity
    }

    static class Range {
        int min
        int max

        Range(int min, int max) {
            this.min = min
            this.max = max
        }

        boolean contains(int value) {
            min <= value && value <= max
        }

        boolean isSingle() {
            max == 1
        }

        @Override
        String toString() {
            "${min}..${max == Integer.MAX_VALUE ? '*' : max}".toString()
        }
    }

}
