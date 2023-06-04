/*
 * Copyright 2013-2023, Seqera Labs
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
 *
 */

package io.seqera.wave.config;

import java.util.List;
import java.util.Map;

/**
 * Spack build options
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
public class SpackOpts {

    final static public String DEFAULT_SPACK_BUILDER_IMAGE = "spack/ubuntu-jammy:v0.19.2";
    final static public String DEFAULT_SPACK_RUNNER_IMAGE = "ubuntu:22.04";
    final static public String DEFAULT_SPACK_OSPACKAGES = "";
    final static public String DEFAULT_SPACK_FLAGS = "-O3";

    public final Boolean checksum;
    public final String builderImage;
    public final String runnerImage;
    public final String osPackages;
    public final String cFlags;
    public final String cxxFlags;
    public final String fFlags;
    public final List<String> commands;

    public SpackOpts() {
        this(Map.of());
    }
    public SpackOpts(Map<String,?> opts) {
        this.checksum = opts.get("checksum") == null || Boolean.parseBoolean(opts.get("checksum").toString());
        this.builderImage = opts.containsKey("builderImage") ? opts.get("builderImage").toString() : DEFAULT_SPACK_BUILDER_IMAGE;
        this.runnerImage = opts.containsKey("runnerImage") ? opts.get("runnerImage").toString() : DEFAULT_SPACK_RUNNER_IMAGE;
        this.osPackages = opts.containsKey("osPackages") ? opts.get("osPackages").toString() : DEFAULT_SPACK_OSPACKAGES;
        this.cFlags = opts.containsKey("cFlags") ? opts.get("cFlags").toString() : DEFAULT_SPACK_FLAGS;
        this.cxxFlags = opts.containsKey("cxxFlags") ? opts.get("cxxFlags").toString() : DEFAULT_SPACK_FLAGS;
        this.fFlags = opts.containsKey("fFlags") ? opts.get("fFlags").toString() : DEFAULT_SPACK_FLAGS;
        this.commands = opts.containsKey("commands") ? (List<String>)opts.get("commands") : null;
    }

}
