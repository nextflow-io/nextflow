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

package nextflow.platform

import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for the per-executor resource label policies
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ResourceLabelPolicyTest extends Specification {

    // -- identity

    def 'should apply the labels verbatim'() {
        given:
        def labels = ['nextflow.io/runName': 'crazy_darwin', 'nextflow.io/repository': 'https://github.com/foo/bar']

        expect:
        ResourceLabelPolicy.IDENTITY.sanitize(labels) == labels
        ResourceLabelPolicy.IDENTITY.sanitizeKey('nextflow.io/runName') == 'nextflow.io/runName'
        ResourceLabelPolicy.IDENTITY.sanitizeValue('https://github.com/foo/bar') == 'https://github.com/foo/bar'
    }

    def 'should handle a null or empty map'() {
        expect:
        POLICY.sanitize(null) == [:]
        POLICY.sanitize([:]) == [:]

        where:
        POLICY << [ResourceLabelPolicy.IDENTITY, ResourceLabelPolicy.AWS, ResourceLabelPolicy.GOOGLE, ResourceLabelPolicy.K8S, ResourceLabelPolicy.AZURE]
    }

    // -- aws

    @Unroll
    def 'should sanitize the aws key'() {
        expect:
        ResourceLabelPolicy.AWS.sanitizeKey(KEY) == EXPECTED

        where:
        KEY                             | EXPECTED
        'nextflow.io/runName'           | 'nextflow.io/runName'
        'seqera.io/platform/workflowId' | 'seqera.io/platform/workflowId'
        'cost centre'                   | 'cost centre'
        'a+b-c=d.e_f:g/h@i'             | 'a+b-c=d.e_f:g/h@i'
        'foo#bar'                       | 'foo_bar'
        'foo(bar)'                      | 'foo_bar_'
        ''                              | ''
        null                            | ''
    }

    def 'should sanitize the aws value'() {
        expect:
        ResourceLabelPolicy.AWS.sanitizeValue('https://github.com/foo/bar') == 'https://github.com/foo/bar'
        ResourceLabelPolicy.AWS.sanitizeValue('some#value') == 'some_value'
    }

    def 'should truncate the aws key and value without a trailing blank'() {
        given:
        def key = 'a' * 127 + ' ' + 'b' * 10
        def value = 'x' * 255 + ' ' + 'y' * 10

        when:
        def k = ResourceLabelPolicy.AWS.sanitizeKey(key)
        def v = ResourceLabelPolicy.AWS.sanitizeValue(value)

        then:
        k == 'a' * 127
        v == 'x' * 255
    }

    def 'should keep an empty aws value'() {
        expect:
        ResourceLabelPolicy.AWS.sanitize(['nextflow.io/revision': '']) == ['nextflow.io/revision': '']
    }

    // -- google

    @Unroll
    def 'should sanitize the google key'() {
        expect:
        ResourceLabelPolicy.GOOGLE.sanitizeKey(KEY) == EXPECTED

        where:
        KEY                             | EXPECTED
        'nextflow.io/runName'           | 'nextflow_io_runname'
        'seqera.io/platform/workflowId' | 'seqera_io_platform_workflowid'
        'foo-bar_baz'                   | 'foo-bar_baz'
        '1foo'                          | 'foo'
        '123'                           | ''
        '-_.'                           | ''
        ''                              | ''
        null                            | ''
    }

    @Unroll
    def 'should sanitize the google value'() {
        expect:
        ResourceLabelPolicy.GOOGLE.sanitizeValue(VALUE) == EXPECTED

        where:
        VALUE                                   | EXPECTED
        'crazy_darwin'                          | 'crazy_darwin'
        'https://github.com/foo/bar'            | 'https___github_com_foo_bar'
        'nf-core/rnaseq'                        | 'nf-core_rnaseq'
        '3.12.0'                                | '3_12_0'
        ''                                      | ''
        null                                    | ''
    }

    def 'should truncate the google key and value'() {
        expect:
        ResourceLabelPolicy.GOOGLE.sanitizeKey('k' * 70) == 'k' * 63
        ResourceLabelPolicy.GOOGLE.sanitizeValue('v' * 70) == 'v' * 63
    }

    def 'should drop the google entries with no key or no value'() {
        expect:
        ResourceLabelPolicy.GOOGLE.sanitize([
            'nextflow.io/runName': 'crazy_darwin',
            '123': 'dropped-no-key',
            'nextflow.io/revision': '' ]) == ['nextflow_io_runname': 'crazy_darwin']
    }

    // -- kubernetes

    @Unroll
    def 'should sanitize the k8s key'() {
        expect:
        ResourceLabelPolicy.K8S.sanitizeKey(KEY) == EXPECTED

        where:
        KEY                             | EXPECTED
        'nextflow.io/runName'           | 'nextflow.io/runName'
        'nextflow.io/resume'            | 'nextflow.io/resume'
        'seqera.io/platform/workflowId' | 'seqera.io/platform_workflowId'
        'runName'                       | 'runName'
        'NEXTFLOW.IO/runName'           | 'nextflow.io/runName'
        'foo bar'                       | 'foo_bar'
        'nextflow.io/-runName-'         | 'nextflow.io/runName'
        'nextflow.io/'                  | ''
        ''                              | ''
        null                            | ''
    }

    @Unroll
    def 'should sanitize the k8s value'() {
        expect:
        ResourceLabelPolicy.K8S.sanitizeValue(VALUE) == EXPECTED

        where:
        VALUE                                       | EXPECTED
        'crazy_darwin'                              | 'crazy_darwin'
        'https://github.com/foo/bar'                | 'github.com_foo_bar'
        'git@github.com:foo/bar.git'                | 'git_github.com_foo_bar.git'
        'ssh://git@github.com/foo/bar'              | 'git_github.com_foo_bar'
        '/leading/and/trailing/'                    | 'leading_and_trailing'
        ''                                          | ''
        null                                        | ''
    }

    def 'should truncate the k8s value leaving an alphanumeric trailing char'() {
        expect:
        ResourceLabelPolicy.K8S.sanitizeValue('a' * 62 + '-' + 'b' * 10) == 'a' * 62
        ResourceLabelPolicy.K8S.sanitizeValue('a' * 70) == 'a' * 63
        ResourceLabelPolicy.K8S.sanitizeKey('nextflow.io/' + 'k' * 70) == 'nextflow.io/' + 'k' * 63
    }

    def 'should drop the k8s entries with no key or no value'() {
        expect:
        ResourceLabelPolicy.K8S.sanitize([
            'nextflow.io/runName': 'crazy_darwin',
            '///': 'dropped-no-key',
            'nextflow.io/repository': '' ]) == ['nextflow.io/runName': 'crazy_darwin']
    }

    // -- azure

    @Unroll
    def 'should sanitize the azure key'() {
        expect:
        ResourceLabelPolicy.AZURE.sanitizeKey(KEY) == EXPECTED

        where:
        KEY                             | EXPECTED
        'nextflow.io/runName'           | 'nextflow.io/runName'
        'seqera.io/platform/workflowId' | 'seqera.io/platform/workflowId'
        'microsoft.foo'                 | 'foo'
        'Microsoft_Bar'                 | 'Bar'
        'microsoft'                     | ''
        ''                              | ''
        null                            | ''
    }

    def 'should keep the azure values as they are'() {
        expect:
        ResourceLabelPolicy.AZURE.sanitize(['nextflow.io/repository': 'https://github.com/foo/bar', 'microsoft': 'dropped', 'nextflow.io/revision': ''])
                == ['nextflow.io/repository': 'https://github.com/foo/bar', 'nextflow.io/revision': '']
    }

    // -- all

    def 'should keep the order of the source map'() {
        given:
        def labels = ['nextflow.io/projectName': 'nf-core/rnaseq', 'nextflow.io/runName': 'crazy_darwin', 'nextflow.io/resume': 'false']

        expect:
        POLICY.sanitize(labels).values() as List == ['nf-core/rnaseq', 'crazy_darwin', 'false'].collect { POLICY.sanitizeValue(it) }

        where:
        POLICY << [ResourceLabelPolicy.IDENTITY, ResourceLabelPolicy.AWS, ResourceLabelPolicy.GOOGLE, ResourceLabelPolicy.K8S, ResourceLabelPolicy.AZURE]
    }

    def 'should have a name'() {
        expect:
        ResourceLabelPolicy.IDENTITY.name == 'identity'
        ResourceLabelPolicy.AWS.name == 'aws'
        ResourceLabelPolicy.GOOGLE.name == 'google'
        ResourceLabelPolicy.K8S.name == 'k8s'
        ResourceLabelPolicy.AZURE.name == 'azure'
    }

}
