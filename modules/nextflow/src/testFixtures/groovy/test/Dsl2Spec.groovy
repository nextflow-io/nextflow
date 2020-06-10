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

package test


import java.nio.file.Path

import groovy.util.logging.Slf4j
import nextflow.NextflowMeta
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Dsl2Spec extends BaseSpec {

    def setupSpec() { NextflowMeta.instance.enableDsl2() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }

    def dsl_eval(String str) {
        new MockScriptRunner().setScript(str).execute()
    }

    def dsl_eval(Path path) {
        new MockScriptRunner().setScript(path).execute()
    }


    def dsl_eval(String entry, String str) {
        new MockScriptRunner()
                .setScript(str).execute(null, entry)
    }
}
