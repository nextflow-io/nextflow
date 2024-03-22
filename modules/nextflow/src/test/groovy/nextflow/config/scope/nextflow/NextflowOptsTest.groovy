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

package nextflow.config.scope.nextflow

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowOptsTest extends Specification {

    def 'should validate equals and hashcode' () {
        given:
        def o1 = new NextflowOpts([defaults: [publishDir: [mode:'foo']]])
        def o2 = new NextflowOpts([defaults: [publishDir: [mode:'foo']]])
        def o3 = new NextflowOpts([defaults: [publishDir: [mode:'bar']]])

        expect:
        o1 == o2
        o1 != o3
        and:
        o1.hashCode() == o2.hashCode()
        o1.hashCode() != o3.hashCode()
    }

    def 'should create empty nextflow opts' () {
        when:
        def nextflow = new NextflowOpts([:])
        then:
        !nextflow.defaults.publishDir.enabled
        !nextflow.defaults.publishDir.mode
        !nextflow.defaults.publishDir.failOnError
        !nextflow.defaults.publishDir.contentType
        !nextflow.defaults.publishDir.overwrite
        !nextflow.defaults.publishDir.contentType
        !nextflow.defaults.publishDir.storageClass
        !nextflow.defaults.publishDir.tags
    }

    def 'should create nextflow publishdir opts' () {
        when:
        def nextflow = new NextflowOpts([
            defaults: [
                publishDir: [
                    mode:'foo',
                    enabled: 'true',
                    failOnError: 'true',
                    contentType: 'some-content',
                    overwrite: 'true',
                    storageClass: 'some-storage',
                    tags: ['this':'one', 'that': 'two']
                ]
            ]])
        then:
        nextflow.defaults.publishDir.mode == 'foo'
        nextflow.defaults.publishDir.enabled
        nextflow.defaults.publishDir.failOnError
        nextflow.defaults.publishDir.contentType == 'some-content'
        nextflow.defaults.publishDir.overwrite
        nextflow.defaults.publishDir.storageClass == 'some-storage'
        nextflow.defaults.publishDir.tags == ['this':'one', 'that': 'two']
    }

}
