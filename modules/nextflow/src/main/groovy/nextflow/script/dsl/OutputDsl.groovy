/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.script.dsl

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.exception.ScriptRuntimeException
import nextflow.script.OutputCollection
import nextflow.script.WorkflowPublisher
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Implements the DSL for top-level workflow outputs
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class OutputDsl {

    private Path path = Path.of('.')

    private List<OutputCollection> collections = []

    void path(String path) {
        this.path = path as Path
    }

    void collect(String name, Closure closure) {
        final dsl = new OutputCollectionDsl()
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()
        this.collections << dsl.build()
    }

    WorkflowPublisher build() {
        new WorkflowPublisher(path, collections)
    }

}

@CompileStatic
class OutputCollectionDsl {

    private String path = '.'

    private List<OutputCollection.Selector> selectors = []

    private OutputCollection.Index index

    void path(String path) {
        this.path = path
    }

    void select(String name) {
        this.selectors << new OutputCollection.Selector(name)
    }

    void select(String name, Closure closure) {
        final dsl = new SelectorDsl()
        dsl.name(name)
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()
        final selector = dsl.build()
        if( selector )
            this.selectors << selector
    }

    void index(Closure closure) {
        final dsl = new IndexDsl()
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()
        this.index = dsl.build()
    }

    OutputCollection build() {
        new OutputCollection(path, selectors, index)
    }

    static class SelectorDsl {
        String name
        boolean enabled = true
        String path = '.'
        String pattern

        void name(String name) {
            this.name = name
        }

        void when(boolean enabled) {
            this.enabled = enabled
        }

        void path(String path) {
            this.path = path
        }

        void pattern(String pattern) {
            this.pattern = pattern
        }

        OutputCollection.Selector build() {
            enabled
                ? new OutputCollection.Selector(name, path, pattern)
                : null
        }
    }

    static class IndexDsl {
        private String format
        private String path

        void format(String format) {
            this.format = format
        }

        void path(String path) {
            this.path = path
        }

        OutputCollection.Index build() {
            new OutputCollection.Index(format, path)
        }
    }

}
