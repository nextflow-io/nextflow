/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
