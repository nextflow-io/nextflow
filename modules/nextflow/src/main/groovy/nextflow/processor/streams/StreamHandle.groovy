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

package nextflow.processor.streams


import groovy.transform.CompileStatic
/**
 * Model a stream handle for process stream input/output definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class StreamHandle {
    enum Direction { IN, OUT }
    Direction direction
    StreamId streamId

    private StreamHandle(Direction dir, StreamId ref) {
        this.direction = dir
        this.streamId = ref
    }

    static StreamHandle input(ref) {
        if( ref instanceof StreamId )
            return new StreamHandle(Direction.IN, ref)
        throw new IllegalArgumentException("Invalid stream refence: $ref")
    }

    static StreamHandle output(ref) {
        if( ref instanceof StreamId )
            return new StreamHandle(Direction.OUT, ref)
        throw new IllegalArgumentException("Invalid stream refence: $ref")
    }

    String name() {
        // This value is replaced into the command script
        // as replacement for the stream variable placeholder
        // it should be a valid linux file name.
        // For this name should be a corresponding named linux named piped
        return "/tmp/nf-${streamId.id()}.stream"
    }

    String toString() {
        return name()
    }
}
