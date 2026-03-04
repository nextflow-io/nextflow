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
package nextflow.script.namespaces;

import java.nio.file.Path;
import java.util.Collection;
import java.util.Map;

import groovy.lang.Closure;
import groovy.transform.NamedParam;
import groovy.transform.NamedParams;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.Namespace;
import nextflow.script.types.Channel;
import nextflow.script.types.Value;

public interface ChannelNamespace extends Namespace {

    @Description("""
        Create a channel that emits nothing.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#empty)
    """)
    Channel<?> empty();

    @Deprecated
    @Description("""
        Create a channel that emits each argument.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#from)
    """)
    <E> Channel<E> from(E... values);

    @Deprecated
    @Description("""
        Create a channel that emits each element in a collection.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#from)
    """)
    <E> Channel<E> from(Collection<E> values);

    @Description("""
        Create a channel that emits all file pairs matching a glob pattern.

        An optional closure can be used to customize the grouping strategy.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
    """)
    Channel<?> fromFilePairs(Map<String,?> opts, String pattern, Closure grouping);

    @Description("""
        Create a channel that emits all paths matching a name or glob pattern.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#frompath)
    """)
    Channel<Path> fromLineage(Map<String,?> opts);

    @Description("""
        Create a channel that emits each element in a collection.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#fromlist)
    """)
    <E> Channel<E> fromList(Collection<E> values);

    @Description("""
        Create a channel that emits all paths matching a name or glob pattern.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#frompath)
    """)
    Channel<Path> fromPath(
        @NamedParams({
            @NamedParam(value = "checkIfExists", type = Boolean.class),
            @NamedParam(value = "followLinks", type = Boolean.class),
            @NamedParam(value = "glob", type = Boolean.class),
            @NamedParam(value = "hidden", type = Boolean.class),
            @NamedParam(value = "maxDepth", type = Integer.class),
            @NamedParam(value = "relative", type = Boolean.class),
            @NamedParam(value = "type", type = String.class),
        })
        Map<String,?> opts,
        String pattern
    );
    Channel<Path> fromPath(String pattern);

    @Description("""
        Create a channel that queries the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) database and emits all FASTQ files matching the given project or accession ids.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#fromsra)
    """)
    Channel<?> fromSRA(Map<String,?> opts, String query);

    @Description("""
        Create a channel that emits an incrementing index (starting from zero) at a periodic interval.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#interval)
    """)
    Channel<Integer> interval(String interval);

    @Description("""
        Create a channel that emits each argument.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#of)
    """)
    <E> Channel<E> of(E... values);

    @Description("""
        Create a channel that emits all values in the given topic.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#topic)
    """)
    Channel<?> topic(String name);

    @Description("""
        Create a value channel.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#value)
    """)
    <V> Value<V> value(V value);

    @Description("""
        Create a channel that watches for filesystem events for all files matching the given pattern.

        [Read more](https://nextflow.io/docs/latest/reference/channel.html#watchpath)
    """)
    Channel<Path> watchPath(String filePattern, String events);

}
