/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.types

import groovy.transform.CompileStatic
import groovy.transform.Immutable

/**
 * Models an instance spot price record
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Immutable
@CompileStatic
class CloudSpotPrice implements Serializable, Cloneable {

    /**
     * The instance type identifier e.g. {@code m4.xlarge}
     */
    String type

    /**
     * The spot price in USD
     */
    String price

    /**
     * The instance availability zone
     */
    String zone

    /**
     * The instance description
     */
    String description

    Date timestamp


}
