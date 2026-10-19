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

package nextflow.file

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

/**
 * The two attributes of a stored object that an {@link ObjectStoreReader} listing carries beside the
 * member's path: its exact byte {@code size} and its last-modified time in {@code mtime} millis.
 *
 * <p>Together they are the change guard of the record-backed file identities: a member whose size and
 * mtime still match the recorded pair is taken to be the same bytes, so its stored token can be
 * reused without reading it. They travel as a value rather than a formatted {@code "size:mtime"}
 * string so that no call site has to parse them back — a parse that could fail, and so needed a
 * sentinel for "unknown" that every caller had to remember to check.
 *
 * <p>{@link #equals} is that guard: the comparison is between whole values, one per member, and no
 * caller reduces them to a string first. {@link #toString} is a debugging aid — it carries no
 * compatibility constraint, and nothing persisted depends on its shape.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@EqualsAndHashCode
class ObjectMeta {

    /** The object's exact size in bytes. */
    final long size

    /** The object's last-modified time, in milliseconds since the epoch. */
    final long mtime

    ObjectMeta(long size, long mtime) {
        this.size = size
        this.mtime = mtime
    }

    /** A compact {@code "size:mtime"} rendering, for logs. Not a wire format — see the class note. */
    @Override
    String toString() {
        return size + ':' + mtime
    }
}
