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
package nextflow.processor.hash

import groovy.transform.CompileStatic
import groovy.transform.Memoized

/**
 * The standard task hash specs, one per historical behaviour of master's hash, loaded
 * from the JSON resources beside this package.
 *
 * These are frozen. Editing a published spec changes the hash of every task already
 * recorded under it; a new behaviour gets a new id and a new file.
 */
@CompileStatic
class StdSpecs {

    static final String RESOURCE_DIR = '/nextflow/processor/hash'

    static final List<String> IDS = ['std/v1', 'std/v2', 'std/v3', 'std/v4']

    static TaskHashSpec getSTD_V1() { byId('std/v1') }
    static TaskHashSpec getSTD_V2() { byId('std/v2') }
    static TaskHashSpec getSTD_V3() { byId('std/v3') }
    static TaskHashSpec getSTD_V4() { byId('std/v4') }

    @Memoized
    static List<TaskHashSpec> all() {
        return IDS.collect { String id -> load(id) }
    }

    static TaskHashSpec byId(String id) {
        final found = all().find { TaskHashSpec it -> it.id == id }
        if( !found ) {
            throw new IllegalArgumentException("Unknown task hash spec: ${id} -- available: ${IDS.join(', ')}")
        }
        return found
    }

    private static TaskHashSpec load(String id) {
        final path = "${RESOURCE_DIR}/${id.replace('/', '-')}.json"
        final stream = StdSpecs.getResourceAsStream(path)
        if( stream == null ) {
            throw new IllegalStateException("Missing task hash spec resource: ${path}")
        }
        try {
            return TaskHashSpecLoader.load(stream.getText('UTF-8'), path)
        }
        finally {
            stream.close()
        }
    }
}
