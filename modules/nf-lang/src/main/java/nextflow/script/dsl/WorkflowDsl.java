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
package nextflow.script.dsl;

import java.util.List;
import java.util.Map;

import groovy.lang.Closure;
import nextflow.script.types.Channel;

/**
 * DSL scope for workflow definitions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public interface WorkflowDsl extends DslScope {

    // support channel operators with pipe (`|`)

    @Operator Object branch(Closure closure);
    @Operator Channel buffer(Closure openingCondition, Closure closingCondition);
    @Operator Channel collate(int size, int step, boolean remainder);
    @Operator Channel collect();
    @Operator Channel collectFile(Map<String,?> opts, Closure closure);
    @Operator Channel combine(Map<String,?> opts, Object right);
    @Operator Channel concat(Channel... others);
    @Operator Channel count();
    @Operator Channel cross(Channel right);
    @Operator Channel distinct();
    // dump not supported by pipe (#4176)
    @Operator Channel filter(Closure<Boolean> condition);
    @Operator Channel first(Object criteria);
    @Operator Channel flatMap(Closure transform);
    @Operator Channel flatten();
    @Operator Channel groupTuple(Map<String,?> opts);
    @Operator Channel ifEmpty(Object value);
    @Operator Channel join(Channel right);
    @Operator Channel last();
    @Operator Channel map(Closure transform);
    @Operator Channel max(Closure comparator);
    @Operator Channel merge(Channel... others);
    @Operator Channel min(Closure comparator);
    @Operator Channel mix(Channel... others);
    @Operator Object multiMap(Closure closure);
    @Operator Channel randomSample(int n, Long seed);
    @Operator Channel reduce(Object seed, Closure accumulator);
    @Operator void set(Closure holder);
    @Operator Channel splitCsv(Map<String,?> opts);
    @Operator Channel splitFasta(Map<String,?> opts);
    @Operator Channel splitFastq(Map<String,?> opts);
    @Operator Channel splitText(Map<String,?> opts, Closure closure);
    @Operator void subscribe(Closure closure);
    @Operator Channel sum(Closure closure);
    @Operator Channel take(int n);
    // tap not supported by pipe (#3970)
    @Operator Channel toList();
    @Operator Channel toSortedList();
    @Operator Channel transpose(Map<String,?> opts);
    @Operator Channel unique(Closure comparator);
    @Operator Channel until(Closure<Boolean> condition);
    @Operator Channel view(Closure transform);

}
