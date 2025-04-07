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
package nextflow.config.dsl;

import java.nio.file.Path;
import java.util.Map;

import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.DslScope;

/**
 * The built-in constants and functions in a config file.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public interface ConfigDsl extends DslScope {

    @Deprecated
    @Constant("baseDir")
    @Description("""
        Alias of `projectDir`.
    """)
    Path getBaseDir();

    @Constant("launchDir")
    @Description("""
        The directory where the workflow was launched.
    """)
    Path getLaunchDir();

    @Constant("params")
    @Description("""
        Map of workflow parameters specified in the config file or as command line options.
    """)
    Map<String,Object> getParams();

    @Constant("projectDir")
    @Description("""
        The directory where the main script is located.
    """)
    Path getProjectDir();

    @Description("""
        Get the value of an environment variable from the launch environment.
    """)
    String env(String name);

}
