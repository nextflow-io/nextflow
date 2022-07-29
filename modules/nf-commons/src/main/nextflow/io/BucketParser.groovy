/*
 * Copyright (c) 2019-2020, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

package nextflow.io

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.StringUtils
/**
 * Parse a cloud bucket uri and decompose in scheme, bucket and path
 * components eg. s3://foo/some/file.txt
 * - scheme: s3
 * - bucket: foo
 * - path  : /some/file.txt
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false, ignoreNulls=false)
@EqualsAndHashCode(includeFields = true)
@CompileStatic
class BucketParser {

    private String scheme
    private String bucket
    private Path path

    String getScheme() { scheme }
    String getBucket() { bucket }
    Path getPath() { path }

    protected BucketParser(String s, String b, String p) {
        this.scheme = s
        this.bucket = b
        this.path = p
                ? (p.startsWith('/') || !bucket ? Path.of(p) : Path.of('/'+p) )
                : Path.of('/')
    }

    protected BucketParser() {
        this.path = Path.of('/')
    }

    static BucketParser from(String uri) {
        new BucketParser().parse(uri)
    }

    BucketParser parse(String uri) {
        final m = StringUtils.URL_PROTOCOL.matcher(uri)
        if( !m.matches() ) {
            path = Path.of(uri)
            return this
        }

        this.scheme = m.group(1)
        final location = m.group(2)

        final p = location.indexOf('/')
        if( p==-1 ) {
            bucket = location
            path = Path.of('/')
        }
        else {
            bucket = location.substring(0,p)
            path = Path.of(location.substring(p))
        }

        if( bucket.startsWith('/') || bucket.endsWith('/') )
            throw new IllegalArgumentException("Invalid bucket URI path: $uri")

        return this
    }


    /**
     * @return The object key, essentially the same as {@link #path} stripping the leading slash character
     */
    String getKey() {
        final result = path.toString()
        return result.substring(1)
    }

}
