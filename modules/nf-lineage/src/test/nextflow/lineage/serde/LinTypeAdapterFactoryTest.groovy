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

package nextflow.lineage.serde


import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.LinModel
import nextflow.lineage.model.v1beta1.TaskRun
import spock.lang.Specification

class LinTypeAdapterFactoryTest extends Specification {

    def 'should read both new and old JSON formats'() {
        given:
        def encoder = new LinEncoder()

        when: 'create old format JSON (without spec wrapper)'
        def oldFormatJson = """
            {
                "version": "${LinModel.VERSION}",
                "kind": "FileOutput",
                "path": "/path/to/file",
                "checksum": {
                    "value": "hash_value",
                    "algorithm": "hash_algorithm",
                    "mode": "standard"
                },
                "source": "lid://source",
                "workflow": "lid://workflow",
                "task": "lid://task",
                "size": 1234
            }
        """

        and: 'deserialize old format'
        def oldResult = encoder.decode(oldFormatJson)

        then: 'should work correctly with backward compatibility'
        oldResult instanceof FileOutput
        oldResult.path == "/path/to/file"
        oldResult.checksum.value == "hash_value"
        oldResult.source == "lid://source"
        oldResult.size == 1234

        when: 'create new format JSON (with spec wrapper)'
        def newFormatJson = """
            {
                "version": "${LinModel.VERSION}",
                "kind": "FileOutput",
                "spec": {
                    "path": "/path/to/file",
                    "checksum": {
                        "value": "hash_value",
                        "algorithm": "hash_algorithm",
                        "mode": "standard"
                    },
                    "source": "lid://source",
                    "workflow": "lid://workflow",
                    "task": "lid://task",
                    "size": 1234
                }
            }
        """

        and: 'deserialize new format'
        def newResult = encoder.decode(newFormatJson)

        then: 'should work correctly'
        newResult instanceof FileOutput
        newResult.path == "/path/to/file"
        newResult.checksum.value == "hash_value"
        newResult.source == "lid://source"
        newResult.size == 1234
    }

    def 'should handle TaskRun old format'() {
        given:
        def encoder = new LinEncoder()

        when: 'create old format TaskRun JSON'
        def oldFormatJson = """
            {
                "version": "${LinModel.VERSION}",
                "kind": "TaskRun",
                "sessionId": "session123",
                "name": "testTask",
                "codeChecksum": {
                    "value": "hash123",
                    "algorithm": "nextflow",
                    "mode": "standard"
                },
                "script": "echo hello",
                "input": [],
                "container": "ubuntu:latest",
                "conda": null,
                "spack": null,
                "architecture": "amd64",
                "globalVars": {},
                "binEntries": []
            }
        """

        and: 'deserialize old format'
        def result = encoder.decode(oldFormatJson)

        then: 'should work correctly'
        result instanceof TaskRun
        result.sessionId == "session123"
        result.name == "testTask"
        result.codeChecksum.value == "hash123"
        result.script == "echo hello"
        result.container == "ubuntu:latest"
        result.architecture == "amd64"
    }

    def 'should reject JSON without version field'() {
        given:
        def encoder = new LinEncoder()

        when: 'try to deserialize JSON without version'
        def invalidJson = """
            {
                "kind": "FileOutput",
                "path": "/path/to/file"
            }
        """
        encoder.decode(invalidJson)

        then: 'should throw exception'
        thrown(Exception)
    }

    def 'should reject JSON with wrong version'() {
        given:
        def encoder = new LinEncoder()

        when: 'try to deserialize JSON with wrong version'
        def invalidJson = """
            {
                "version": "wrong/version",
                "kind": "FileOutput",
                "path": "/path/to/file"
            }
        """
        encoder.decode(invalidJson)

        then: 'should throw exception'
        thrown(Exception)
    }

    def 'should reject JSON without kind field'() {
        given:
        def encoder = new LinEncoder()

        when: 'try to deserialize JSON without kind'
        def invalidJson = """
            {
                "version": "${LinModel.VERSION}",
                "path": "/path/to/file"
            }
        """
        encoder.decode(invalidJson)

        then: 'should throw exception'
        thrown(Exception)
    }
}
