package io.nextflow.gradle.util

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * S3 bucket tokenizer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode(includeFields = true)
@CompileStatic
class BucketTokenizer {

    static final public Pattern URL_PROTOCOL = ~/^([a-zA-Z0-9]*):\\/\\/(.+)/

    private String scheme
    private String bucket
    private String path
    private boolean directory

    String getScheme() { scheme }
    String getBucket() { bucket }
    String getPath() { path }
    boolean isDirectory() { directory }

    protected BucketTokenizer(String s, String b, String p, boolean dir=false) {
        this.scheme = s
        this.bucket = b
        this.path = p
        this.directory = dir
    }

    protected BucketTokenizer() {}

    static BucketTokenizer from(String uri) {
        new BucketTokenizer().parse(uri)
    }

    BucketTokenizer parse(String uri) {
        final m = URL_PROTOCOL.matcher(uri)
        if( !m.matches() ) {
            return this
        }

        this.scheme = m.group(1)
        final location = m.group(2)

        final p = location.indexOf('/')
        if( p==-1 ) {
            bucket = location
            path = ''
        }
        else {
            bucket = location.substring(0,p)
            path = location.substring(p)
        }

        directory = path.endsWith('/')

        if( bucket.startsWith('/') || bucket.endsWith('/') )
            throw new IllegalArgumentException("Invalid bucket URI path: $uri")

        while(path.endsWith('/'))
            path = path.substring(0,path.length()-1)

        return this
    }

    BucketTokenizer withPath(String newPath) {
        final dir = newPath?.endsWith('/')
        while(newPath.endsWith('/'))
            newPath = newPath.substring(0,newPath.length()-1)
        new BucketTokenizer(this.scheme, this.bucket, newPath, dir)
    }

    String toString() {
        def result = scheme
        if( bucket )
            result += '://' + bucket
        if( bucket && path )
            result += '/' + path
        return result
    }

    /**
     * @return The object key, essentially the same as {@link #path} stripping the leading slash character
     */
    String getKey() {
        path?.startsWith('/') ? path.substring(1) : path
    }

}
