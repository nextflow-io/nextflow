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
 */

package test

import spock.lang.Specification

import groovy.util.logging.Slf4j

/**
 * Base specification class - It wraps each test into begin-close "test name" strings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BaseSpec extends Specification {

    def setup() {
        log.info "TEST BEGIN [${specificationContext.currentIteration.name}]"
    }

    def cleanup() {
        log.info "TEST CLOSE [${specificationContext.currentIteration.name}]"
    }


}
