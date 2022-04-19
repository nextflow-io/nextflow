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

package com.upplication.s3fs.experiment

import com.upplication.s3fs.experiment.AtomicBigInteger
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AtomicBigIntegerTest extends Specification {

    def 'should increment and get' ( ) {
        given:
        def acc = new AtomicBigInteger()

        expect:
        and:
        acc.incrementAndGet() == 1
        acc.incrementAndGet() == 2
        and:
        acc.getAndIncrement() == 2
        acc.getAndIncrement() == 3
        and:
        acc.get() == 4

    }

    def 'should increment one and get' ( ) {
        given:
        def acc = new AtomicBigInteger()

        expect:
        and:
        acc.getAndIncrement(1) == 0
        and:
        acc.get() == 1

    }

}
