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

package nextflow.cloud.google.batch

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpServer
import io.seqera.cloudinfo.client.CloudInfoClient
import nextflow.cloud.google.batch.GoogleBatchMachineTypeSelector.MachineType
import nextflow.cloud.types.PriceModel
import nextflow.util.MemoryUnit
import spock.lang.Specification

class GoogleBatchMachineTypeSelectorTest extends Specification {

    static final MACHINE_TYPES = [
            new MachineType(type: 'e2-type01', family: 'e2', 'spotPrice': 0.001, 'onDemandPrice': 0.01, 'cpusPerVm': 1, 'memPerVm': 1),
            new MachineType(type: 'n1-type02', family: 'n1', 'spotPrice': 0.002, 'onDemandPrice': 0.02, 'cpusPerVm': 2, 'memPerVm': 2),
            new MachineType(type: 'e2-type03', family: 'e2', 'spotPrice': 0.010, 'onDemandPrice': 0.05, 'cpusPerVm': 4, 'memPerVm': 4),
            new MachineType(type: 'n2-type04', family: 'n2', 'spotPrice': 0.011, 'onDemandPrice': 0.15, 'cpusPerVm': 4, 'memPerVm': 4),
            new MachineType(type: 'e2-type05', family: 'e2', 'spotPrice': 0.020, 'onDemandPrice': 0.20, 'cpusPerVm': 6, 'memPerVm': 6),
            new MachineType(type: 'n1-type06', family: 'n1', 'spotPrice': 0.025, 'onDemandPrice': 0.25, 'cpusPerVm': 6, 'memPerVm': 6),
            new MachineType(type: 'm1-type07', family: 'm1', 'spotPrice': 0.030, 'onDemandPrice': 0.30, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'm2-type08', family: 'm2', 'spotPrice': 0.036, 'onDemandPrice': 0.35, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'n2-type09', family: 'n2', 'spotPrice': 0.040, 'onDemandPrice': 0.40, 'cpusPerVm': 10, 'memPerVm': 10),
            new MachineType(type: 'c2-type10', family: 'c2', 'spotPrice': 0.045, 'onDemandPrice': 0.45, 'cpusPerVm': 10, 'memPerVm': 10),
            new MachineType(type: 'c4-type11', family: 'c4', 'spotPrice': 0.040, 'onDemandPrice': 0.40, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'c4a-type12', family: 'c4a', 'spotPrice': 0.038, 'onDemandPrice': 0.38, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'c4d-type13', family: 'c4d', 'spotPrice': 0.039, 'onDemandPrice': 0.39, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'n4-type14', family: 'n4', 'spotPrice': 0.035, 'onDemandPrice': 0.35, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'n4a-type15', family: 'n4a', 'spotPrice': 0.033, 'onDemandPrice': 0.33, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'n4d-type16', family: 'n4d', 'spotPrice': 0.034, 'onDemandPrice': 0.34, 'cpusPerVm': 8, 'memPerVm': 8),
    ]

    def 'should select best machine type'() {
        given:
        final selector = Spy(GoogleBatchMachineTypeSelector) {
            getAvailableMachineTypes(REGION, SPOT) >> MACHINE_TYPES
        }
        expect:
        selector.bestMachineType(CPUS, MEM, REGION, SPOT, FUSION, FAMILIES).type == EXPECTED

        where:
        CPUS | MEM  | REGION | SPOT  | FUSION | FAMILIES                   | EXPECTED
        1    | 1000 | 'reg'  | true  | false  | null                       | 'e2-type01'
        1    | 1000 | 'reg'  | false | true   | null                       | 'n1-type02'
        4    | 4000 | 'reg'  | false | false  | []                         | 'e2-type03'
        4    | 4000 | 'reg'  | true  | false  | []                         | 'n2-type04'
        4    | 4000 | 'reg'  | false | true   | []                         | 'n2-type04'
        6    | 6000 | 'reg'  | true  | false  | null                       | 'e2-type05'
        6    | 6000 | 'reg'  | true  | true   | null                       | 'n1-type06'
        6    | 6000 | 'reg'  | true  | false  | ['n1-*', 'm1-*']           | 'n1-type06'
        8    | 8000 | 'reg'  | true  | false  | null                       | 'm1-type07'
        8    | 8000 | 'reg'  | false | false  | ['m?-*', 'c2-*']           | 'm2-type08'
        8    | 8000 | 'reg'  | false | false  | ['m1-type07', 'm2-type66'] | 'm1-type07'
        8    | 8000 | 'reg'  | true  | false  | ['c4-*']                   | 'c4-type11'
        8    | 8000 | 'reg'  | true  | false  | ['c4a-*']                  | 'c4a-type12'
        8    | 8000 | 'reg'  | true  | false  | ['c4d-*']                  | 'c4d-type13'
        8    | 8000 | 'reg'  | true  | false  | ['n4-*']                   | 'n4-type14'
        8    | 8000 | 'reg'  | true  | false  | ['n4a-*']                  | 'n4a-type15'
        8    | 8000 | 'reg'  | true  | false  | ['n4d-*']                  | 'n4d-type16'


    }

    def 'should not select a machine type'() {
        given:
        final selector = Spy(GoogleBatchMachineTypeSelector) {
            getAvailableMachineTypes(REGION, SPOT) >> MACHINE_TYPES
        }
        when:
        selector.bestMachineType(CPUS, MEM, REGION, SPOT, SSD, FAMILIES)

        then:
        thrown(NoSuchElementException)

        where:
        CPUS | MEM  | REGION | SPOT | SSD   | FAMILIES
        8    | 9000 | 'reg'  | true | true  | ['s2-*', 'e2-*']
        12   | 1000 | 'reg'  | true | false | null

    }

    // a trimmed response of the Seqera cloud info products API for google/europe-west2
    static final String CLOUD_INFO_PRODUCTS = '''
        {"products": [
          {"category": "General purpose", "type": "e2-small", "family": "e2", "onDemandPrice": 0.06373944000000001,
           "spotPrice": [{"zone": "europe-west2-a", "price": 0.025}, {"zone": "europe-west2-b", "price": 0.026}, {"zone": "europe-west2-c", "price": 0.027}],
           "cpusPerVm": 2, "memPerVm": 2, "gpusPerVm": 0, "ntwPerf": "1 Gbit/s", "ntwPerfCategory": "low",
           "zones": ["europe-west2-c", "europe-west2-b", "europe-west2-a"],
           "attributes": {"cpu": "2", "instanceTypeCategory": "General purpose", "memory": "2", "networkPerfCategory": "low"},
           "currentGen": false, "features": ["burst", "sched"], "maximumPersistentDisks": 16, "maximumPersistentDisksSizeGb": 3072,
           "isSharedCpu": true, "bootDiskTypes": ["pd-balanced", "pd-ssd", "pd-standard"]},
          {"category": "General purpose", "type": "n2-standard-4", "family": "n2", "onDemandPrice": 0.250248,
           "spotPrice": [{"zone": "europe-west2-a", "price": 0.052256000000000004}, {"zone": "europe-west2-b", "price": 0.052256000000000004}],
           "cpusPerVm": 4, "memPerVm": 16, "gpusPerVm": 0, "ntwPerf": "10 Gbit/s", "ntwPerfCategory": "high",
           "zones": ["europe-west2-c", "europe-west2-b", "europe-west2-a"],
           "attributes": {"cpu": "4", "instanceTypeCategory": "General purpose", "memory": "16", "networkPerfCategory": "high"},
           "currentGen": false, "features": ["sched"], "attachableLocalSsdCounts": [1, 2, 4, 8, 16, 24],
           "maximumPersistentDisks": 128, "maximumPersistentDisksSizeGb": 263168, "bootDiskTypes": ["pd-balanced", "pd-ssd", "pd-standard"]},
          {"category": "General purpose", "type": "a3-edgegpu-8g", "family": "a3", "onDemandPrice": 100.91807079189041,
           "spotPrice": [],
           "cpusPerVm": 208, "memPerVm": 1872, "gpusPerVm": 8, "ntwPerf": "400 Gbit/s", "ntwPerfCategory": "extra",
           "zones": ["europe-west2-c", "europe-west2-b", "europe-west2-a"],
           "attributes": {"cpu": "208", "instanceTypeCategory": "General purpose", "memory": "1872", "networkPerfCategory": "extra"},
           "currentGen": false, "features": ["gpu", "h100-80gb", "nvidia", "sched", "ssd", "x86"],
           "gpuModel": "nvidia-h100", "acceleratorType": "nvidia-h100-80gb", "localSsdPartitions": 16, "localSsdGiB": 6000}
        ],
        "scrapingTime": 1794411332182}
        '''

    def 'should parse Seqera cloud info API'() {
        given: 'a local server replying with a recorded cloud info response'
        def requestPath = null
        def server = HttpServer.create(new InetSocketAddress(InetAddress.loopbackAddress, 0), 0)
        server.createContext('/') { HttpExchange exchange ->
            requestPath = exchange.requestURI.path
            final body = CLOUD_INFO_PRODUCTS.getBytes('UTF-8')
            exchange.responseHeaders.add('Content-Type', 'application/json')
            exchange.sendResponseHeaders(200, body.length)
            exchange.responseBody.withCloseable { it.write(body) }
        }
        server.start()
        and:
        def selector = new GoogleBatchMachineTypeSelector()
        selector.@cloudInfo = CloudInfoClient.builder()
            .endpoint("http://localhost:${server.address.port}")
            .build()

        when:
        def types = selector.getAvailableMachineTypes('europe-west2', true)

        then:
        requestPath == '/api/v1/providers/google/services/compute/regions/europe-west2/products'
        types.size() == 3
        and:
        with(types[0]) {
            type == 'e2-small'
            family == 'e2'
            location == 'europe-west2'
            priceModel == PriceModel.spot
            Math.abs(spotPrice - 0.026f) < 1e-6
            Math.abs(onDemandPrice - 0.06373944f) < 1e-6
            cpusPerVm == 2
            memPerVm == 2
            gpusPerVm == 0
        }
        and:
        types[1].type == 'n2-standard-4'
        types[1].family == 'n2'
        Math.abs(types[1].spotPrice - 0.052256f) < 1e-6
        types[1].memPerVm == 16
        and:
        types[2].type == 'a3-edgegpu-8g'
        types[2].family == 'a3'
        types[2].spotPrice == 0.0f
        types[2].cpusPerVm == 208
        types[2].gpusPerVm == 8

        cleanup:
        server?.stop(0)
    }

    def 'should find first valid disk size'() {
        expect:
        GoogleBatchMachineTypeSelector.INSTANCE.findFirstValidSize(MemoryUnit.of(REQUESTED), ALLOWED) == MemoryUnit.of(EXPECTED)

        where:
        REQUESTED | ALLOWED   | EXPECTED
        '100 GB'  | [2, 4, 8] | '750 GB'
        '100 GB'  | [1, 2, 4] | '375 GB'
        '500 GB'  | [1, 2]    | '750 GB'
        '1 TB'    | [1, 2]    | '750 GB'
    }

    def 'should find valid local disk size given the machine type'() {
        expect:
        final machineType = new MachineType(type: TYPE, family: FAMILY, cpusPerVm: CPUS)
        GoogleBatchMachineTypeSelector.INSTANCE.findValidLocalSSDSize(MemoryUnit.of(REQUESTED), machineType) == MemoryUnit.of(EXPECTED)

        where:
        REQUESTED | TYPE              | FAMILY | CPUS | EXPECTED
        '100 GB'  | 'n1-highmem-8'    | 'n1'   | 8    | '375 GB'
        '375 GB'  | 'n2-highcpu-16'   | 'n2'   | 16   | '750 GB'
        '780 GB'  | 'n2d-standard-48' | 'n2d'  | 48   | '1500 GB'
        '200 GB'  | 'c2-standard-4'   | 'c2'   | 4    | '375 GB'
        '50 GB'   | 'c2d-highmem-56'  | 'c2d'  | 56   | '1500 GB'
        '750 GB'  | 'm3-megamem-64'   | 'm3'   | 64   | '1500 GB'
        '100 GB'  | 'c4-standard-8-lssd'  | 'c4'   | 8    | '0'
        '100 GB'  | 'c4a-standard-8-lssd' | 'c4a'  | 8    | '0'
        '100 GB'  | 'c4d-standard-8-lssd' | 'c4d'  | 8    | '0'
    }

    def 'should know when hyperdisk is required'() {
        expect:
        GoogleBatchMachineTypeSelector.INSTANCE.isHyperdiskOnly(TYPE) == EXPECTED

        where:
        TYPE                | EXPECTED
        'c4-standard-8'     | true
        'c4a-standard-8'    | true
        'c4d-standard-8'    | true
        'n4-standard-8'     | true
        'n4a-standard-8'    | true
        'n4d-standard-8'    | true
        'n1-standard-8'     | false
        'n2-standard-8'     | false
        'e2-standard-8'     | false
        'c2-standard-8'     | false
    }

    def 'should know when to install GPU drivers'() {
        expect:
        final machineType = new MachineType(type: TYPE, gpusPerVm: GPUS)
        GoogleBatchMachineTypeSelector.INSTANCE.installGpuDrivers(machineType) == EXPECTED

        where:
        TYPE            | GPUS | EXPECTED
        'n2-standard-4' | 0    | false
        'n2-standard-4' | 1    | true
        'a2-highgpu-1g' | 0    | true
        'a3-highgpu-1g' | 0    | true
        'g2-standard-4' | 0    | true
    }

    def 'should detect non-configurable local SSD'() {
        expect:
        final machineType = new MachineType(type: TYPE, family: FAMILY)
        GoogleBatchMachineTypeSelector.INSTANCE.notConfigurableLocalSSD(machineType) == EXPECTED

        where:
        TYPE                        | FAMILY | EXPECTED
        // c3/c3d with -lssd suffix → true
        'c3-standard-8-lssd'        | 'c3'   | true
        'c3d-standard-8-lssd'       | 'c3d'  | true
        // c4/c4a/c4d with -lssd suffix → true
        'c4-standard-8-lssd'        | 'c4'   | true
        'c4a-standard-8-lssd'       | 'c4a'  | true
        'c4d-standard-8-lssd'       | 'c4d'  | true
        // a3 family → always true regardless of type
        'a3-highgpu-8g'             | 'a3'   | true
        'a3-megagpu-64g'            | 'a3'   | true
        // a2-ultragpu- prefix → true regardless of family
        'a2-ultragpu-1g'            | 'a2'   | true
        'a2-ultragpu-8g'            | 'a2'   | true
        // c3/c4 without -lssd suffix → false
        'c3-standard-8'             | 'c3'   | false
        'c4-standard-8'             | 'c4'   | false
        // a2 non-ultragpu → false
        'a2-highgpu-1g'             | 'a2'   | false
        // unrelated families → false
        'n2-standard-4'             | 'n2'   | false
        'e2-standard-8'             | 'e2'   | false
    }
}
