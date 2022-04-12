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

import com.upplication.s3fs.experiment.TransferRateMeter
import org.junit.Ignore
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Ignore
class TransferRateMeterTest extends Specification {

    def 'should compute rate' () {
        given:
        def meter = new TransferRateMeter()

        when:
        def t = Thread.start { 100.times { index -> sleep 100; meter.inc(10_000) } }
        t.join()
        then:
        noExceptionThrown()

    }
}
