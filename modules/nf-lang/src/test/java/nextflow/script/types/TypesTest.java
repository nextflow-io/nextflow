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
package nextflow.script.types;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import nextflow.script.dsl.Types;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.GenericsType;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.MethodSource;

import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies how types are rendered and compared.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypesTest {

    @Test
    void renderType() {
        assertThat(Types.getName(ClassHelper.MAP_TYPE)).isEqualTo("Map<K, V>");

        var cn = new ClassNode(Map.class);
        assertThat(Types.getName(cn)).isEqualTo("Map");

        cn.setGenericsTypes(new GenericsType[] {
            new GenericsType(ClassHelper.STRING_TYPE),
            new GenericsType(ClassHelper.STRING_TYPE)
        });
        assertThat(Types.getName(cn)).isEqualTo("Map<String, String>");
    }

    // resolves to nextflow.script.types.Record, which shadows java.lang.Record
    static class Sample implements Record {}

    static Stream<Arguments> runtimeTypes() {
        return Stream.of(
            arguments(Path.of("/a").getClass(), "Path"),
            arguments(ArrayList.class, "List"),
            arguments(Sample.class, "Record"),
            arguments(Record.class, "Record"),
            arguments(String.class, "String")
        );
    }

    @ParameterizedTest
    @MethodSource("runtimeTypes")
    void renderRuntimeType(Class<?> type, String name) {
        assertThat(Types.getName(type)).isEqualTo(name);
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', textBlock = """
        getExt   | Map<String, String>
        getLabel | List<String>
        """)
    void renderMethodReturnType(String method, String name) {
        var cn = ClassHelper.makeCached(TaskConfig.class);
        var returnType = cn.getDeclaredMethods(method).get(0).getReturnType();
        assertThat(Types.getName(returnType)).isEqualTo(name);
    }

    @Test
    void renderFunctionalType() {
        var cn = ClassHelper.makeCached(Channel.class);

        var filter = cn.getDeclaredMethods("filter").get(0).getParameters()[0];
        assertThat(Types.getName(filter.getType())).isEqualTo("(E) -> Boolean");

        var map = cn.getDeclaredMethods("map").get(0).getParameters()[0];
        assertThat(Types.getName(map.getType())).isEqualTo("(E) -> R");

        // reduce is overloaded and declared-method order is unspecified
        var reduceParams = cn.getDeclaredMethods("reduce").get(0).getParameters();
        var reduce = reduceParams[reduceParams.length - 1];
        assertThat(Types.getName(reduce.getType())).isIn("(E, E) -> E", "(R, E) -> R");
    }

    static Stream<Arguments> assignableTypes() {
        return Stream.of(
            arguments(Integer.class, Integer.class, true),
            arguments(Integer.class, Float.class, false),
            arguments(Float.class, Integer.class, true),
            arguments(Float.class, Float.class, true),

            arguments(String.class, String.class, true),
            arguments(String.class, Integer.class, false),
            arguments(Integer.class, String.class, false),

            arguments(Iterable.class, Bag.class, true),
            arguments(Iterable.class, List.class, true),
            arguments(Iterable.class, Set.class, true)
        );
    }

    @ParameterizedTest
    @MethodSource("assignableTypes")
    void determineWhetherTypeIsAssignableToAnother(Class<?> target, Class<?> source, boolean result) {
        assertThat(Types.isAssignableFrom(ClassHelper.makeCached(target), ClassHelper.makeCached(source))).isEqualTo(result);
    }

}
