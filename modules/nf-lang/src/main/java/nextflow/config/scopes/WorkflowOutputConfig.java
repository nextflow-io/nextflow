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
package nextflow.config.scopes;

import java.util.Map;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class WorkflowOutputConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        *Currently only supported for S3.*

        Specify the media type a.k.a. [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_Types) of published files (default: `false`). Can be a string (e.g. `'text/html'`), or `true` to infer the content type from the file extension.
    """)
    public Object contentType;

    @ConfigOption
    @Description("""
        Enable or disable publishing (default: `true`).
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        When `true`, the workflow will not fail if a file can't be published for some reason (default: `false`).
    """)
    public boolean ignoreErrors;

    @ConfigOption
    @Description("""
        The file publishing method (default: `'symlink'`).
    """)
    public String mode;

    @ConfigOption
    @Description("""
        When `true` any existing file in the specified folder will be overwritten (default: `'standard'`).
    """)
    public Object overwrite;

    @ConfigOption
    @Description("""
        *Currently only supported for S3.*

        Specify the storage class for published files.
    """)
    public String storageClass;

    @ConfigOption
    @Description("""
        *Currently only supported for S3.*

        Specify arbitrary tags for published files.
    """)
    public Map<String,String> tags;

}
