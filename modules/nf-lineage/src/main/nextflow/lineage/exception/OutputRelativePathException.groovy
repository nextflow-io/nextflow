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
package nextflow.lineage.exception

/**
 * Exception to indicate the an output path is not relative to the output dir.
 * It is used to detect the cases where publishDir is used with Data Lineage.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class OutputRelativePathException extends Exception {
}
