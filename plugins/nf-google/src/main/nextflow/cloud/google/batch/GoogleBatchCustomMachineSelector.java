package nextflow.cloud.google.batch;

import com.google.common.collect.ImmutableMap;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

/**
 * Base on Google Cloud documentation example
 * https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type#java
 */
class GoogleBatchCustomMachineSelector {

    // This class defines the configurable parameters for a custom VM.
    static final class TypeLimits {

        int[] allowedCores;
        int minMemPerCore;
        int maxMemPerCore;
        int extraMemoryLimit;
        boolean allowExtraMemory;

        TypeLimits(int[] allowedCores, int minMemPerCore, int maxMemPerCore, boolean allowExtraMemory,
                   int extraMemoryLimit) {
            this.allowedCores = allowedCores;
            this.minMemPerCore = minMemPerCore;
            this.maxMemPerCore = maxMemPerCore;
            this.allowExtraMemory = allowExtraMemory;
            this.extraMemoryLimit = extraMemoryLimit;
        }
    }

    enum CpuSeries {
        N1("custom"),
        N2("n2-custom"),
        N2D("n2d-custom"),
        E2("e2-custom"),
        E2_MICRO("e2-custom-micro"),
        E2_SMALL("e2-custom-small"),
        E2_MEDIUM("e2-custom-medium");

        private static final Map<String, CpuSeries> ENUM_MAP;

        static {
            ENUM_MAP = init();
        }

        // Build an immutable map of String name to enum pairs.
        static Map<String, CpuSeries> init() {
            Map<String, CpuSeries> map = new ConcurrentHashMap<>();
            for (CpuSeries instance : CpuSeries.values()) {
                map.put(instance.getCpuSeries(), instance);
            }
            return Collections.unmodifiableMap(map);
        }

        private final String cpuSeries;

        CpuSeries(String cpuSeries) {
            this.cpuSeries = cpuSeries;
        }

        static CpuSeries get(String name) {
            return ENUM_MAP.get(name);
        }

        String getCpuSeries() {
            return this.cpuSeries;
        }
    }

    // Returns the array of integers within the given range, incremented by the specified step.
    // start (inclusive): starting number of the range
    // stop (inclusive): ending number of the range
    // step : increment value
    static int[] getNumsInRangeWithStep(int start, int stop, int step) {
        return IntStream.range(start, stop).filter(x -> (x - start) % step == 0).toArray();
    }

    static int gbToMb(int value) {
        return value << 10;
    }

    static int[] concat(int[] a, int[] b) {
        int[] result = new int[a.length + b.length];
        System.arraycopy(a, 0, result, 0, a.length);
        System.arraycopy(b, 0, result, a.length, b.length);
        return result;
    }

    // This enum correlates a machine type with its limits.
    // The limits for various CPU types are described in:
    // https://cloud.google.com/compute/docs/general-purpose-machines
    enum Limits {
        CPUSeries_E2(new TypeLimits(getNumsInRangeWithStep(2, 33, 2), 512, 8192, false, 0)),
        CPUSeries_E2MICRO(new TypeLimits(new int[]{}, 1024, 2048, false, 0)),
        CPUSeries_E2SMALL(new TypeLimits(new int[]{}, 2048, 4096, false, 0)),
        CPUSeries_E2MEDIUM(new TypeLimits(new int[]{}, 4096, 8192, false, 0)),
        CPUSeries_N2(
                new TypeLimits(concat(getNumsInRangeWithStep(2, 33, 2), getNumsInRangeWithStep(36, 129, 4)),
                        512, 8192, true, gbToMb(624))),
        CPUSeries_N2D(
                new TypeLimits(new int[]{2, 4, 8, 16, 32, 48, 64, 80, 96}, 512, 8192, true, gbToMb(768))),
        CPUSeries_N1(
                new TypeLimits(concat(new int[]{1}, getNumsInRangeWithStep(2, 97, 2)), 922, 6656, true,
                        gbToMb(624)));

        private final TypeLimits typeLimits;

        Limits(TypeLimits typeLimits) {
            this.typeLimits = typeLimits;
        }

        TypeLimits getTypeLimits() {
            return typeLimits;
        }
    }

    static ImmutableMap<String, Limits> typeLimitsMap = ImmutableMap.<String, Limits>builder()
            .put("N1", Limits.CPUSeries_N1)
            .put("N2", Limits.CPUSeries_N2)
            .put("N2D", Limits.CPUSeries_N2D)
            .put("E2", Limits.CPUSeries_E2)
            .put("E2_MICRO", Limits.CPUSeries_E2MICRO)
            .put("E2_SMALL", Limits.CPUSeries_E2SMALL)
            .put("E2_MEDIUM", Limits.CPUSeries_E2SMALL)
            .build();


    static String customMachineTypeUri(String cpuSeries, int coreCount, int memory) {

        if (!Arrays.asList(CpuSeries.E2.cpuSeries, CpuSeries.N1.cpuSeries, CpuSeries.N2.cpuSeries,
                CpuSeries.N2D.cpuSeries).contains(cpuSeries)) {
            throw new Error(String.format("Incorrect cpu type: %s", cpuSeries));
        }

        TypeLimits typeLimit = Objects.requireNonNull(typeLimitsMap.get(CpuSeries.get(cpuSeries).name())).typeLimits;

        // Perform the following checks to verify if the requested parameters are allowed.
        // Find more information about limitations of custom machine types at:
        // https://cloud.google.com/compute/docs/general-purpose-machines#custom_machine_types

        // 1. Check the number of cores and if the coreCount is present in allowedCores.
        final int cores = coreCount;
        if (typeLimit.allowedCores.length > 0 && Arrays.stream(typeLimit.allowedCores)
                .noneMatch(x -> x == cores)) {

            // Find nearest
            coreCount = findValidCoreCount(cpuSeries, typeLimit, cores);
        }

        // 2. Memory must be a multiple of 256 MB
        if (memory % 256 != 0) {
            int m = memory / 256;
            memory = (m + 1) * 256;
        }

        // 3. Check if the requested memory isn't too little
        if (memory < coreCount * typeLimit.minMemPerCore) {
            // Add more memory to match the minimum limit
            memory = coreCount * typeLimit.minMemPerCore;
        }

        // 4. Check if the requested memory isn't too much
        if (memory > coreCount * typeLimit.maxMemPerCore && !typeLimit.allowExtraMemory) {
            // Add more cores if needed
            coreCount = findValidCoreCount(cpuSeries, typeLimit, memory / typeLimit.maxMemPerCore);
        }

        // 5. Check if the requested memory isn't too large
        if (memory > typeLimit.extraMemoryLimit && typeLimit.allowExtraMemory) {
            // Set it to the maximum
            memory = typeLimit.extraMemoryLimit;
            //TODO fail or warn ??
        }

        // Check if the CPU Series is E2 and return the custom machine type in the form of a string
        // acceptable by Compute Engine API.
        if (Arrays.asList(CpuSeries.E2_SMALL.cpuSeries, CpuSeries.E2_MICRO.cpuSeries,
                CpuSeries.E2_MEDIUM.cpuSeries).contains(cpuSeries)) {
            return String.format("%s-%s", cpuSeries, memory);
        }

        // Check if extended memory was requested and return the extended custom machine type
        // in the form of a string acceptable by Compute Engine API.
        if (memory > coreCount * typeLimit.maxMemPerCore) {
            return String.format("%s-%s-%s-ext", cpuSeries, coreCount, memory);
        }

        // Return the custom machine type in the form of a standard string
        // acceptable by Compute Engine API.
        return String.format("%s-%s-%s", cpuSeries, coreCount, memory);
    }

    private static int findValidCoreCount(String cpuSeries, TypeLimits typeLimit, int cores) {
        return Arrays.stream(typeLimit.allowedCores).boxed()
                .min(Comparator.comparingInt(i -> Math.abs(i - cores)))
                .orElseThrow(() -> new Error(String.format(
                        "Invalid number of cores requested. "
                                + "Number of cores requested for CPU %s should be one of: %s",
                        cpuSeries,
                        Arrays.toString(typeLimit.allowedCores))));
    }
}
