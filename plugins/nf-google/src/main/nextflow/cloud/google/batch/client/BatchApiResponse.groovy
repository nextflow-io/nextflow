/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.client

import groovy.transform.CompileStatic


/**
 * Model a Batch API response
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchApiResponse {

    private int code

    private InputStream stream

    private String text

    BatchApiResponse(int code, InputStream stream) {
        this.code = code
        this.stream = stream
    }

    String toString() {
        "code=$code; stream=$stream"
    }

    int getCode() { code }

    InputStream getStream() { stream }

    String getText() {
        if( text == null ) {
            text = stream?.text
        }
        return text
    }
}
