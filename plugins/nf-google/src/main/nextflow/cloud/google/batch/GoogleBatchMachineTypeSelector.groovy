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

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Immutable
import groovy.transform.Memoized

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
    private static final List<String> DEFAULT_FAMILIES_WITH_SSD = ['n1-*', 'n2-*', 'n2d-*', 'c2-*', 'c2d-*', 'm3-*']

    @Immutable
    static class MachineType {
        String type
        String family
        float spotPrice
        float onDemandPrice
        int cpusPerVm
        int memPerVm
    }

    String bestMachineType(int cpus, int memoryMB, String region, boolean spot, boolean localSSD, List<String> families) {
        final machineTypes = getAvailableMachineTypes(region)
        if (families == null)
            families = []

        // Check if a specific machine type was defined
        if (families.size() == 1) {
            final familyOrType = families.get(0)
            if (familyOrType.contains("custom-"))
                return familyOrType

            if (machineTypes.find { it.type == familyOrType })
                return familyOrType
        }

        final memoryGB = Math.ceil(memoryMB / 1024.0 as float) as int

        // Use only families that can have a local SSD
        if (!families && localSSD)
            families = DEFAULT_FAMILIES_WITH_SSD

        // All types are valid if no families are defined, otherwise at least it has to start with one of the given values
        final matchMachineType = {String t -> !families || families.find { matchType(it, t) }}

        // find machines with enough resources and SSD local disk
        final validMachineTypes = machineTypes.findAll {
                    it.cpusPerVm >= cpus &&
                    it.memPerVm >= memoryGB &&
                    matchMachineType(it.type)
        }.collect()

        final sortedByCost = validMachineTypes.sort {
            (it.cpusPerVm > 2 || it.memPerVm > 2 ? FAMILY_COST_CORRECTION.get(it.family, 1.0) : 1.0) * (spot ? it.spotPrice : it.onDemandPrice)
        }

        return sortedByCost.first().type
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
    protected List<MachineType> getAvailableMachineTypes(String region) {
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
                    memPerVm: it.memPerVm as int
            )
        }
    }
}
