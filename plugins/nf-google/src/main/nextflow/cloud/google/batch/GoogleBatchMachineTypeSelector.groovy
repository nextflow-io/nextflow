/*
 * Copyright 2023, Seqera Labs
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
 *
 */
package nextflow.cloud.google.batch

import java.math.RoundingMode

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Immutable
import groovy.transform.Memoized
import nextflow.cloud.types.PriceModel
import nextflow.util.MemoryUnit

/**
 * Choose best machine type that fits the requested resources and
 * reduces the estimated cost on current location.
 *
 * Average spot and on-demand prices are requested once per execution
 * from Seqera cloud info API.
 *
 * If cloud info service is not available it will fallback to use
 * the user provided machine type or default to Google Batch automatic
 * selection from the requested resource.
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
@CompileStatic
class GoogleBatchMachineTypeSelector {

    static GoogleBatchMachineTypeSelector INSTANCE = new GoogleBatchMachineTypeSelector()

    private static final CLOUD_INFO_API = "https://cloudinfo.seqera.io/api/v1"

    /*
     * Some families CPUs are faster so this is a cost correction factor
     * for processes that request more than 2 CPUs or 2GB, smaller processes
     * we assume that do not have high CPU usage.
     * https://cloud.google.com/compute/docs/cpu-platforms
     */
    private static final Map<String, BigDecimal> FAMILY_COST_CORRECTION = [
            'e2' : 1.0,   // Mix of processors, tend to be similar in performance to N1

            // INTEL
            'n1' : 1.0,   // Skylake, Broadwell, Haswell, Sandy Bridge, and Ivy Bridge ~2.7 Ghz
            'n2' : 0.85,  // Intel Xeon Gold 6268CL ~3.4 Ghz
            'c2' : 0.8,   // Intel Xeon Gold 6253CL ~3.8 Ghz
            'm1' : 1.0,   // Intel Xeon E7-8880V4 ~2.7 Ghz
            'm2' : 0.85,  // Intel Xeon Platinum 8280L ~3.4 Ghz
            'm3' : 0.85,  // Intel Xeon Platinum 8373C ~3.4 Ghz
            'a2' : 0.9,   // Intel Xeon Platinum 8273CL ~2.9 Ghz

            // AMD
            't2d': 1.0,   // AMD EPYC Milan ~2.7 Ghz
            'n2d': 1.0,   // AMD EPYC Milan ~2.7 Ghz
    ]

    /*
     * Families that will be use as default if Fusion is enabled but no list is provided
     * https://cloud.google.com/compute/docs/disks#local_ssd_machine_type_restrictions
     * LAST UPDATE 2023-02-17
     */
    private static final List<String> DEFAULT_FAMILIES_FOR_FUSION = ['n1-*', 'n2-*', 'n2d-*', 'c2-*', 'c2d-*', 'm3-*']

    private static final List<String> DEFAULT_FAMILIES = ['n1-*', 'n2-*', 'n2d-*', 'c2-*', 'c2d-*', 'm1-*', 'm2-*', 'm3-*', 'e2-*']

    @Immutable
    static class MachineType {
        String type
        String family
        String location
        float spotPrice
        float onDemandPrice
        int cpusPerVm
        int memPerVm
        PriceModel priceModel
    }

    MachineType bestMachineType(int cpus, int memoryMB, String region, boolean spot, boolean fusionEnabled, List<String> families) {
        final machineTypes = getAvailableMachineTypes(region, spot)
        if (families == null)
            families = Collections.<String>emptyList()

        // Check if a specific machine type was defined
        if (families.size() == 1) {
            final familyOrType = families.get(0)
            if (familyOrType.contains("custom-"))
                return new MachineType(type: familyOrType, family: 'custom', cpusPerVm: cpus, memPerVm: memoryMB, location: region, priceModel: spot ? PriceModel.spot : PriceModel.standard)

            final machineType = machineTypes.find { it.type == familyOrType }
            if( machineType )
                return machineType
        }

        final memoryGB = Math.ceil(memoryMB / 1024.0 as float) as int

        if (!families ) {
            families = fusionEnabled
                    ? DEFAULT_FAMILIES_FOR_FUSION
                    : DEFAULT_FAMILIES
        }

        // All types are valid if no families are defined, otherwise at least it has to start with one of the given values
        final matchMachineType = {String type -> !families || families.find { matchType(it, type) }}

        // find machines with enough resources and SSD local disk
        final validMachineTypes = machineTypes.findAll {
                    it.cpusPerVm >= cpus &&
                    it.memPerVm >= memoryGB &&
                    matchMachineType(it.type)
        }.collect()

        final sortedByCost = validMachineTypes.sort {
            (it.cpusPerVm > 2 || it.memPerVm > 2 ? FAMILY_COST_CORRECTION.get(it.family, 1.0) : 1.0) * (spot ? it.spotPrice : it.onDemandPrice)
        }

        return sortedByCost.first()
    }

    protected boolean matchType(String family, String vmType) {
        if (!family)
            return true
        if (family.contains('*'))
            family = family.toLowerCase().replaceAll(/\*/, '.*')
        if (family.contains('?'))
            family = family.toLowerCase().replaceAll(/\?/, '.{1}')

        return vmType =~ /(?i)^${family}$/
    }

    @Memoized
    protected List<MachineType> getAvailableMachineTypes(String region, boolean spot) {
        final priceModel = spot ? PriceModel.spot : PriceModel.standard
        final json = "${CLOUD_INFO_API}/providers/google/services/compute/regions/${region}/products".toURL().text
        final data = new JsonSlurper().parseText(json)
        final products = data['products'] as List<Map>
        final averageSpotPrice = (List<Map> prices) -> prices.collect{it.price as float}.average() as float

        products.collect {
            new MachineType(
                    type: it.type,
                    family: it.type.toString().split('-')[0],
                    spotPrice: averageSpotPrice(it.spotPrice as List<Map>),
                    onDemandPrice: it.onDemandPrice as float,
                    cpusPerVm: it.cpusPerVm as int,
                    memPerVm: it.memPerVm as int,
                    location: region,
                    priceModel: priceModel
            )
        }
    }

    /**
     * Find valid local SSD size. See: https://cloud.google.com/compute/docs/disks#local_ssd_machine_type_restrictions
     *
     * @param requested Amount of disk requested
     * @param machineType Machine type
     * @return Next greater multiple of 375 GB that is a valid size for the given machine type
     */
    protected MemoryUnit findValidLocalSSDSize(MemoryUnit requested, MachineType machineType) {

        if( machineType.family == "n1" )
            return findFirstValidSize(requested, [1,2,3,4,5,6,7,8,16,24])

        if( machineType.family == "n2" ) {
            if( machineType.cpusPerVm < 12 )
                return findFirstValidSize(requested, [1,2,4,8,16,24])
            if( machineType.cpusPerVm < 22 )
                return findFirstValidSize(requested, [2,4,8,16,24])
            if( machineType.cpusPerVm < 42 )
                return findFirstValidSize(requested, [4,8,16,24])
            if( machineType.cpusPerVm < 82 )
                return findFirstValidSize(requested, [8,16,24])
            return findFirstValidSize(requested, [16,24])
        }

        if( machineType.family == "n2d" ) {
            if( machineType.cpusPerVm < 32 )
                return findFirstValidSize(requested, [1,2,4,8,16,24])
            if( machineType.cpusPerVm < 64 )
                return findFirstValidSize(requested, [2,4,8,16,24])
            if( machineType.cpusPerVm < 96 )
                return findFirstValidSize(requested, [4,8,16,24])
            return findFirstValidSize(requested, [8,16,24])
        }

        if( machineType.family == "c2" ) {
            if( machineType.cpusPerVm < 16 )
                return findFirstValidSize(requested, [1,2,4,8])
            if( machineType.cpusPerVm < 30 )
                return findFirstValidSize(requested, [2,4,8])
            if( machineType.cpusPerVm < 60 )
                return findFirstValidSize(requested, [4,8])
            return findFirstValidSize(requested, [8])
        }

        if( machineType.family == "c2d" ) {
            if( machineType.cpusPerVm < 32 )
                return findFirstValidSize(requested, [1,2,4,8])
            if( machineType.cpusPerVm < 56 )
                return findFirstValidSize(requested, [2,4,8])
            if( machineType.cpusPerVm < 112 )
                return findFirstValidSize(requested, [4,8])
            return findFirstValidSize(requested, [8])
        }

        if( machineType.family == "m3" ) {
            if ( machineType.type == 'm3-megamem-128' || machineType.type == 'm3-ultramem-128' )
                return findFirstValidSize(requested, [8])
            return findFirstValidSize(requested, [4,8])
        }

        // other special families the user must provide a valid size
        return requested
    }

    /**
     * Find first valid disk size given the possible mounted partition
     *
     * @param requested Requested disk size
     * @param allowedPartitions Valid number of disks of 375.GB.
     * @return
     */
    protected MemoryUnit findFirstValidSize(MemoryUnit requested, List<Integer> allowedPartitions) {

        // Sort the possible number of disks
        allowedPartitions.sort()

        // Minimum number of 375.GB disks to fulfill the requested size
        final disks = (requested.toGiga() / 375).setScale(0, RoundingMode.UP).toInteger()

        // Find first valid number of disks
        def numberOfDisks = allowedPartitions.find { it >= disks}
        if( !numberOfDisks )
            numberOfDisks = allowedPartitions.last()

        return new MemoryUnit( numberOfDisks * 375L * (1<<30) )
    }

}
