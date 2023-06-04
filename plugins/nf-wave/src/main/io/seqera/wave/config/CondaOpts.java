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
 * Conda build options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CondaOpts {
    final public static String DEFAULT_MAMBA_IMAGE = "mambaorg/micromamba:1.4.2";

    final public String mambaImage;
    final public List<String> commands;
    final public String basePackages;

    public CondaOpts(Map<String,?> opts) {
        this.mambaImage = opts.containsKey("mambaImage") ? opts.get("mambaImage").toString(): DEFAULT_MAMBA_IMAGE;
        this.commands = opts.containsKey("commands") ? (List<String>)opts.get("commands") : null;
        this.basePackages = (String)opts.get("basePackages");
    }

}
