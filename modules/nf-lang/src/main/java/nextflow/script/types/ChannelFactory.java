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
package nextflow.script.types;

import java.nio.file.Path;
import java.util.Collection;
import java.util.Map;

import groovy.lang.Closure;

public interface ChannelFactory {

    Channel empty();

    <E> Channel<E> from(E... values);

    <E> Channel<E> from(Collection<E> values);

    Channel fromFilePairs(Map<String,?> opts, String pattern, Closure grouping);

    <E> Channel<E> fromList(Collection<E> values);

    Channel<Path> fromPath(Map<String,?> opts, String pattern);

    Channel fromSRA(Map<String,?> opts, String query);

    <E> Channel<E> of(E... values);

    Channel topic(String name);

    <E> Channel<E> value(E value);

    Channel<Path> watchPath(String filePattern, String events);

}
