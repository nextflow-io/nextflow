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

package nextflow.serde


import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.serde.gson.GsonEncoder
import nextflow.serde.gson.RuntimeTypeAdapterFactory
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MyEncoder extends GsonEncoder<JsonSerializable> {

    MyEncoder() {
        withTypeAdapterFactory(
            RuntimeTypeAdapterFactory.of(JsonSerializable.class, "@type")
                .registerSubtype(Dog.class, "Dog")
                .registerSubtype(Cat.class, "Cat")
        )
    }

}

@EqualsAndHashCode
class Dog implements JsonSerializable {
    private final String name;
    int barkVolume;

    Dog(String name, int barkVolume) {
        this.name = name;
        this.barkVolume = barkVolume;
    }

    String getName() {
        return name;
    }
}

@EqualsAndHashCode
class Cat implements JsonSerializable {
    private final String name;
    boolean likesSun;

    Cat(String name, boolean likesSun) {
        this.name = name;
        this.likesSun = likesSun;
    }

    String getName() {
        return name;
    }
}
