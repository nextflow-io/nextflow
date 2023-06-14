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

    /**
     * Custom Dockerfile `RUN` commands that can be used to customise the target container
     */
    public final List<String> commands;

    /**
     * Spack packages that should be added to any Spack environment requested via Wave
     */
    public final String basePackages;

    public SpackOpts() {
        this(Map.of());
    }
    public SpackOpts(Map<String,?> opts) {
        this.commands = opts.containsKey("commands") ? (List<String>)opts.get("commands") : null;
        this.basePackages = opts.containsKey("basePackages") ? opts.get("basePackages").toString() : null;
    }

}
