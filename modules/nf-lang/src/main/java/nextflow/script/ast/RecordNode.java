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
package nextflow.script.ast;

import java.lang.reflect.Modifier;

import nextflow.script.types.Record;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;

/**
 * A record type definition.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class RecordNode extends ClassNode {
    public RecordNode(String name) {
        super(name, Modifier.PUBLIC | Modifier.FINAL, ClassHelper.OBJECT_TYPE);
        setInterfaces(new ClassNode[] { ClassHelper.makeCached(Record.class) });
    }

    @Override
    public Class getTypeClass() {
        return Record.class;
    }
}
