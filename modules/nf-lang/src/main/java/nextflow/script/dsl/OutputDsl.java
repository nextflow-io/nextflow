/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.dsl;

import java.util.List;
import java.util.Map;

import groovy.lang.Closure;

/**
 * DSL scope for the workflow output definition.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public interface OutputDsl extends DslScope {

    @Constant("args")
    @Description("""
        List of positional arguments specified on the command line.
    """)
    List<String> getArgs();

    @Constant("params")
    @Description("""
        Map of workflow parameters specified in the config file or as command line options.
    """)
    Map<String,Object> getParams();

    @Description("""
        *Currently only supported for S3.*

        Specify the media type a.k.a. [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_Types) of published files (default: `false`). Can be a string (e.g. `'text/html'`), or `true` to infer the content type from the file extension.
    """)
    /* String | Boolean */
    void contentType(Object value);

    @Description("""
        Enable or disable publishing (default: `true`).
    """)
    void enabled(boolean value);

    @Description("""
        When `true`, the workflow will not fail if a file can't be published for some reason (default: `false`).
    """)
    void ignoreErrors(boolean value);

    @Description("""
        Create an index file of the values that were published.
    """)
    void index(Closure closure);

    @Description("""
        The file publishing method (default: `'symlink'`).
    """)
    void mode(String value);

    @Description("""
        When `true` any existing file in the specified folder will be overwritten (default: `'standard'`).
    """)
    /* String | Boolean */
    void overwrite(Object value);

    @Description("""
        Specify the publish path relative to the output directory (default: the target name).
    """)
    void path(String value);

    @Description("""
        *Currently only supported for S3.*

        Specify the storage class for published files.
    """)
    void storageClass(String value);

    @Description("""
        *Currently only supported for S3.*

        Specify arbitrary tags for published files.
    """)
    void tags(Map<String,String> value);

    interface IndexDsl extends DslScope {

        @Description("""
            When `true`, the keys of the first record are used as the column names (default: `false`). Can also be a list of column names.
        """)
        /* List<String> | Boolean */
        void header(Object value);

        @Description("""
            Closure which defines how to transform each published value into a CSV record. The closure should return a list or map. By default, no transformation is applied.
        """)
        void mapper(Closure value);

        @Description("""
            The name of the index file relative to the target path (required).
        """)
        void path(String value);

        @Description("""
            The character used to separate values (default: `','`).
        """)
        void sep(String value);

    }

}
