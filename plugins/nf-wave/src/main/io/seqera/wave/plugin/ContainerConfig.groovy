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

package io.seqera.wave.plugin

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.util.CacheHelper

/**
 * Model a container configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
@ToString(includePackage = false, includeNames = true, ignoreNulls = true)
class ContainerConfig {

    /**
     * Prefixes of env variables configuring the Fusion driver operation (logging,
     * tracing, snapshots) that do not change what a task computes, therefore
     * excluded from the task fingerprint. See {@link #taskHashKey()}
     */
    private static final List<String> SKIP_HASHING_ENV = List.of(
        'FUSION_SHOW_GLOBAL_INFO',
        'FUSION_TRACER_',
        'FUSION_SNAPSHOT_',
        'FUSION_LOG_' )

    List<String> entrypoint
    List<String> cmd
    List<String> env
    String workingDir

    List<ContainerLayer> layers

    ContainerConfig appendLayer(ContainerLayer it)  {
        if( layers==null && it!=null )
            layers = new ArrayList<>(10)
        layers.add(it)
        return this
    }

    ContainerConfig prependLayer(ContainerLayer it)  {
        if( layers==null && it!=null )
            layers = new ArrayList<>(10)
        layers.add(0, it)
        return this
    }

    ContainerConfig plus( ContainerConfig that ) {
        new ContainerConfig(
                entrypoint: that.entrypoint ?: this.entrypoint,
                cmd: that.cmd ?: this.cmd,
                env: mergeEnv(this.env, that.env),
                workingDir: that.workingDir ?: this.workingDir,
                layers: mergeLayers(this.layers, that.layers) )
    }

    protected List<String> mergeEnv(List<String> left, List<String> right) {
        if( left==null && right==null )
            return null
        final result = new ArrayList<String>(10)
        if( left )
            result.addAll(left)
        if( !right )
            return result
        // add the 'right' env to the result
        for(String it : right) {
            final pair = it.tokenize('=')
            // remove existing var because the right list has priority
            final p = result.findIndexOf { it.startsWith(pair[0]+'=')}
            if( p!=-1 ) {
                result.remove(p)
                result.add(p, it)
            }
            else
                result.add(it)
        }
        return result
    }

    protected List<ContainerLayer> mergeLayers(List<ContainerLayer> left, List<ContainerLayer> right) {
        if( left==null && right==null )
            return null
        final result = new ArrayList<ContainerLayer>()
        if( left ) result.addAll(left)
        if( right ) result.addAll(right)
        return result
    }

    String fingerprint() {
        return fingerprint0(env)
    }

    /**
     * Same as {@link #fingerprint()} ignoring the env variables in {@link #SKIP_HASHING_ENV}.
     * Used as the container term of the task hash.
     */
    String taskHashKey() {
        return fingerprint0(env?.findAll { String it -> !SKIP_HASHING_ENV.any { String p -> it.startsWith(p) } })
    }

    private String fingerprint0(List<String> env) {
        final allMeta = new ArrayList()
        allMeta.add( entrypoint ?: 'no-entry' )
        allMeta.add( cmd ?: 'no-cmd' )
        allMeta.add( env ?: 'no-env' )
        allMeta.add( workingDir ?: 'no-workdir')
        final layers0 = layers ?: Collections.<ContainerLayer>emptyList()

        for( ContainerLayer it : layers0 ) {
            if( !it.skipHashing )
                allMeta.add(it.fingerprint())
        }
        return CacheHelper.hasher(allMeta).hash().toString()
    }

    boolean asBoolean() {
        return !empty()
    }

    boolean empty() {
        return !entrypoint &&
            !cmd &&
            !env &&
            !workingDir &&
            !layers
    }
}
