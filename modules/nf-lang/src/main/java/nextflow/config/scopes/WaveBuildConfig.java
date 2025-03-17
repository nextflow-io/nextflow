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

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class WaveBuildConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The container repository where images built by Wave are uploaded.
    """)
    public String repository;

    @ConfigOption
    @Description("""
        The container repository used to cache image layers built by the Wave service.
    """)
    public String cacheRepository;

    public WaveBuildCondaConfig conda;

    public WaveBuildSpackConfig spack;

}

class WaveBuildCondaConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        One or more Conda packages to be always added in the resulting container (default: `conda-forge::procps-ng`).
    """)
    public String basePackages;

    @ConfigOption
    @Description("""
        One or more commands to be added to the Dockerfile used to build a Conda based image.
    """)
    public String commands;

    @ConfigOption
    @Description("""
        The Mamba container image that is used to build Conda based container.
    """)
    public String mambaImage;

}

class WaveBuildSpackConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        One or more Spack packages to be always added in the resulting container.
    """)
    public String basePackages;

    @ConfigOption
    @Description("""
        One or more commands to be added to the Dockerfile used to build a Spack based image.
    """)
    public String commands;

}
