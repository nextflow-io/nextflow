/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.model.v1beta1

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.util.CacheHelper
/**
 * Models a checksum including the value as well as the algortihm and mode used to compute it.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io
 */
@Canonical
@CompileStatic
class Checksum {
    String value
    String algorithm
    String mode

    static Checksum of(String value, String algorithm, CacheHelper.HashMode mode) {
        new Checksum(value, algorithm, mode.toString().toLowerCase())
    }

    static Checksum ofNextflow(String value) {
        new Checksum(CacheHelper.hasher(value).hash().toString(), 'nextflow', CacheHelper.HashMode.DEFAULT().toString().toLowerCase())
    }

    static Checksum ofNextflow(Path path) {
        new Checksum(CacheHelper.hasher(path).hash().toString(), 'nextflow', CacheHelper.HashMode.DEFAULT().toString().toLowerCase())
    }
}
