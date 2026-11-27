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

package io.seqera.config

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * An entry of the {@code seqera.executor.additionalRegions} setting: a region the scheduler
 * can use besides the primary one, optionally with the network to run tasks in.
 *
 * An entry is either a region name (e.g. {@code 'us-east-1'}) or a map such as
 * {@code [region: 'us-east-1', network: [vpcId: 'vpc-1', subnetIds: ['subnet-1'], securityGroups: ['sg-1']]]}.
 */
@ToString(includeNames = true, ignoreNulls = true)
@EqualsAndHashCode
@CompileStatic
class AdditionalRegionOpts {

    static final String SETTING = 'seqera.executor.additionalRegions'

    private static final Set<String> REGION_KEYS = Set.of('region', 'network')

    private static final Set<String> NETWORK_KEYS = Set.of('vpcId', 'subnetIds', 'securityGroups')

    final String region

    /** The VPC ID, or {@code null} when not specified */
    final String vpcId

    /** The subnet IDs, or {@code null} when not specified */
    final List<String> subnetIds

    /** The security group IDs, or {@code null} when not specified */
    final List<String> securityGroups

    AdditionalRegionOpts(String region, String vpcId = null, List<String> subnetIds = null, List<String> securityGroups = null) {
        this.region = region
        this.vpcId = vpcId
        this.subnetIds = subnetIds
        this.securityGroups = securityGroups
    }

    /**
     * Whether a network was specified for this region.
     */
    boolean hasNetwork() {
        return vpcId || subnetIds || securityGroups
    }

    /**
     * Parse the {@code additionalRegions} setting.
     *
     * @param value a comma-separated string, or a list whose entries are region names or maps
     * @return the parsed regions, or {@code null} when none is specified
     */
    static List<AdditionalRegionOpts> parse(Object value) {
        if( value == null )
            return null
        final List entries = value instanceof Collection
            ? new ArrayList(value as Collection)
            : value.toString().tokenize(',') as List
        final result = new ArrayList<AdditionalRegionOpts>()
        for( Object entry : entries ) {
            final opts = parseEntry(entry)
            if( opts )
                result.add(opts)
        }
        return result ?: null
    }

    private static AdditionalRegionOpts parseEntry(Object entry) {
        if( entry == null )
            return null
        if( entry instanceof Map )
            return fromMap(entry as Map)
        final region = entry.toString().trim()
        return region ? new AdditionalRegionOpts(region) : null
    }

    private static AdditionalRegionOpts fromMap(Map map) {
        checkKeys(map, REGION_KEYS, SETTING)
        final region = map.region?.toString()?.trim()
        if( !region )
            throw new IllegalArgumentException("Missing 'region' in '${SETTING}' entry: ${map}")
        final network = map.network
        if( network == null )
            return new AdditionalRegionOpts(region)
        if( network !instanceof Map )
            throw new IllegalArgumentException("Invalid 'network' in '${SETTING}' entry for region '${region}' - expected a map, got: ${network}")
        final net = network as Map
        checkKeys(net, NETWORK_KEYS, "${SETTING}.network")
        return new AdditionalRegionOpts(
            region,
            net.vpcId?.toString()?.trim() ?: null,
            toList(net.subnetIds),
            toList(net.securityGroups) )
    }

    private static void checkKeys(Map map, Set<String> allowed, String setting) {
        for( Object key : map.keySet() ) {
            if( key?.toString() !in allowed )
                throw new IllegalArgumentException("Unknown option '${key}' in '${setting}' - valid options: ${allowed.toSorted().join(', ')}")
        }
    }

    private static List<String> toList(Object value) {
        if( value == null )
            return null
        final List items = value instanceof Collection
            ? new ArrayList(value as Collection)
            : value.toString().tokenize(',') as List
        final result = new ArrayList<String>()
        for( Object it : items ) {
            final str = it?.toString()?.trim()
            if( str )
                result.add(str)
        }
        return result ?: null
    }
}
