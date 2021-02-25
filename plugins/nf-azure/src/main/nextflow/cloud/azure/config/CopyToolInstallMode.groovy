package nextflow.cloud.azure.config

/**
 * Define the strategy to install the azcopy tool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
enum CopyToolInstallMode {
    node,
    task
}
