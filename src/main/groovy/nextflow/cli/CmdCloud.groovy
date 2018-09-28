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

package nextflow.cli
import static nextflow.Const.ROLE_MASTER
import static nextflow.Const.ROLE_WORKER
import static nextflow.cli.CmdHelper.fixEqualsOp
import static nextflow.cloud.CloudConst.TAG_CLUSTER_NAME

import java.nio.file.Paths

import ch.grengine.Grengine
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.CloudConfig
import nextflow.cloud.CloudDriver
import nextflow.cloud.CloudDriverFactory
import nextflow.cloud.types.CloudInstance
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudSpotPrice
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.ui.TableBuilder
import nextflow.ui.TextLabel
import nextflow.util.SysHelper
/**
 * Implements the `cloud` command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
//@CompileStatic -- Don't use. It causes a weird exception: CmdCloud$LaunchMaster cannot be cast to nextflow.cli.CmdCloud
@Parameters(commandDescription = "Manage Nextflow clusters in the cloud")
class CmdCloud extends CmdBase implements UsageAware {

    interface SubCmd {
        String getName()
        String getShortcut()
        void apply(String arg)
        void usage(List<String> result)
    }

    static public final String NAME = 'cloud'

    private List<SubCmd> commands = []

    @Parameter(names=['-driver'], description = 'The name of a cloud driver')
    String driverName

    @Parameter(names=['-t','-instance-type'], description = 'Instance type')
    String instanceType

    @Parameter(names=['-i','-image-id'], description = 'Image ID')
    String imageId

    @Parameter(names=['-c', '-instance-count'], description = 'Instances count')
    int instanceCount

    @Parameter(names='-spot-price', description = 'Price for spot/preemptive instances')
    String spotPrice

    @Parameter(names=['-p','-profile'], description='Configuration profile')
    String profile

    @Parameter(names=['-sort'], description = 'Sort price history by the specified field: type, price, zone, description')
    String sort

    @Parameter(names=['-F','-filter'], description = 'Filter price history by the specified field: type, price, zone, description')
    String filter

    @Parameter(names='-history', description = 'Print price history for the specified instance type')
    String history

    @Parameter(names='-all', description = 'Print all prices and availability zones', arity = 0)
    boolean all

    @Parameter(names=['-r','-region'], description = 'The region to use. Overrides config/env settings.')
    String region

    @Parameter(names='-y', description = 'Skip launch confirmation')
    boolean yes

    @Parameter
    List<String> args

    private CloudDriver driver

    private Map config

    CmdCloud() {
        commands << new CreateCluster()
        commands << new JoinCluster()
        commands << new ListNodes()
        commands << new SpotPrices()
        commands << new Shutdown()
    }

    /**
     * @return The command name i.e. `cloud`
     */
    @Override
    String getName() {
        return NAME
    }

    /**
     * Initialise the command required data structures i.e.
     * 1) Load the cloud driver
     * 2) Read the nextflow config file containing the cloud configuration
     * 3) Check the command specified is valid
     *
     */
    private void init() {
        final availDrivers = CloudDriverFactory.getDriverNames()
        if( !availDrivers )
            throw new AbortOperationException("No cloud drivers are available")

        // no driver was specified -- choose the first available
        if( !driverName ) {
            if( availDrivers.size()>1 ) throw new AbortOperationException("No cloud driver was specified -- Use option -driver to choose one of the following: ${availDrivers.collect { ',' }}")
            driverName = availDrivers[0]
        }
        // check that an existing driver was specified
        else if( !availDrivers.contains(driverName) ) {
            def matches = availDrivers.closest(driverName)
            def msg = "Unknow cloud driver name: $driverName"
            if( matches )
                msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
            throw new AbortOperationException(msg)
        }

        // -- create the config object
        this.config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(Paths.get('.'))
                .setProfile(profile)
                .build()

        Global.setConfig(config)
        this.driver = CloudDriverFactory.getDriver(driverName, [region: this.region])
    }

    /**
     * Command entry point
     */
    @Override
    void run() {
        // -- if no args are provided just show the usage and exit
        if( !args ) {
            usage()
            return
        }

        // -- initialize
        init()

        // -- find out the sub-command
        final cmd = findCommand(args)

        // invoke it
        cmd.apply( args.size()>1 ? args[1] : null)

    }

    /**
     * Creates the {@link CloudConfig} object with the given command line options
     *
     * @param clusterName The name to be assigned to the cloud cluster
     * @return An {@link CloudConfig} object instance
     */
    protected CloudConfig makeConfig(String clusterName) {
        def result = CloudConfig.create(config)
        result.setClusterName(clusterName)
        if( driverName ) result.setDriverName(driverName)
        if( instanceType) result.setInstanceType(instanceType)
        if( imageId ) result.setImageId(imageId)
        if( spotPrice ) result.setSpotPrice(spotPrice)

        result.build()
        result.validate(driver)
        return result
    }

    /**
     * Prints the cloud configuration in a human readable manner
     *
     * @param cloudConfig A {@link CloudConfig} object
     * @param count the number of instances to launch
     */
    protected void printConfig(CloudConfig cloudConfig, int count) {
        println "> cluster name: $cloudConfig.clusterName"
        println "> instances count: $count"
        println "> Launch configuration:"
        println cloudConfig.prettyPrint().indent(  )
    }

    /**
     * Prints the list of worker instances
     *
     * @param instanceIds The list of instance IDs to be shown
     */
    protected void printWorkerInstances( List<String> instanceIds ) {

        // -- show the available nodes
        driver.eachInstanceWithIds(instanceIds) { item ->
            println "  ${item.id}\t  ${item.address}"
        }

        println ''
    }

    /**
     * Prints the information about the cluster master instance
     *
     * @param instanceId The ID of the master instance
     * @param cloudConfig The {@link CloudConfig} object of the current cluster
     */
    protected void printMasterInstance( String instanceId, CloudConfig cloudConfig ) {
        // -- notify that's ready
        driver.eachInstanceWithIds([instanceId]) { CloudInstance it ->
            println "Login in the master node using the following command: \n  ${getSshCommand(cloudConfig, it.address)}"
        }
        println ""
    }

    /**
     * A string describing the SSH command to login into the cluster master instance
     *
     * @param cfg The {@link CloudConfig} cluster configuration
     * @param address The cluster master name public address (public DNS name)
     * @return The SSH command string
     */
    protected String getSshCommand(CloudConfig cfg, String address) {
        if( !address )
            throw new IllegalStateException("Missing instance public DNS name")

        String key
        if( cfg.keyName ) {
            key = "<path to `${cfg.keyName}` key file>"
        }
        else if( cfg.privateKeyFile?.exists() ) {
            key = cfg.privateKeyFile
        }
        else {
            key = "<path to your private key file>"
        }
       "ssh -i ${key} ${cfg.userName}@${address}"
    }


    /**
     * Add one or more node to the specified cluster
     *
     * @param clusterName
     * @param count
     * @return
     */
    private List<String> addNodes0(int count, CloudConfig cloudConfig) {
        assert count, 'Parameter count must be greater than zero'

        // -- validate before launch nodes
        if( !cloudConfig.role )
            throw new IllegalStateException("Missing cloud config `role` attribute")

        // -- submit the launch request
        final nodes = count==1 ? 'node' : 'nodes'
        print "Launching ${cloudConfig.role} $nodes -- Waiting for `running` status.. "
        final instanceIds = driver.launchInstances(count, cloudConfig)

        // -- wait for running status
        driver.waitInstanceStatus(instanceIds, CloudInstanceStatus.STARTED)

        // -- tagging
        driver.tagInstances(instanceIds, cloudConfig)

        if( cloudConfig.role == ROLE_MASTER ) {
            driver.waitInstanceStatus(instanceIds, CloudInstanceStatus.READY)
        }

        println "ready."

        return instanceIds
    }

    /**
     * Prompt the user for confirmation
     *
     * @param message The message to show when asking the user confirmation
     * @throws AbortOperationException when the user answer NO
     */
    private void prompt(String message) {
        print "$message [y/n] "
        char answer = System.in.read()
        if( answer != 'y' && answer != 'Y' ) throw new AbortOperationException()
    }

    /**
     * Validate a cluster name.
     *
     * @param name A string representing the cluster name provided by the user
     */
    @PackageScope
    void checkName(String name) {
        if( !name )
            throw new AbortOperationException("Missing cluster name -- See `nextflow cloud -h` for usage")

        if( name.isNumber() )
            throw new AbortOperationException("Specified cluster name is not valid: $name -- It must be an alphanumeric string")

        def matches = name ==~ /[a-zA-Z][a-zA-Z0-9\-_]+/
        if( !matches ) {
            throw new AbortOperationException("Specified cluster name is not valid: `$name` -- It can only contain alphanumeric, - and _ characters")
        }
    }

    @PackageScope
    Set<String> getClusterNames() {
        def tags = [:]; tags[TAG_CLUSTER_NAME] = null
        def names = new HashSet()
        driver.eachInstanceWithTags(tags) { CloudInstance it -> names.add(it.clusterName) }
        return names
    }

    private SubCmd findCommand(List<String> args) {

        def cmd = commands.find { it.name == args[0] || it.shortcut == args[0] || (it.name.endsWith('s') && it.name[0..-2]==args[0]) }
        if( cmd ) {
            return cmd
        }

        def matches = commands.collect{ it.name }.closest(args[0])
        def msg = "Unknown cloud sub-command: ${args[0]}"
        if( matches )
            msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
        throw new AbortOperationException(msg)
    }

    /**
     * Print the command usage help
     */
    void usage() {
        usage(args)
    }

    /**
     * Print the command usage help
     *
     * @param args The arguments as entered by the user
     */
    void usage(List<String> args) {

        def result = []
        if( !args ) {
            result << this.getClass().getAnnotation(Parameters).commandDescription()
            result << 'Usage: nextflow cloud <command> [arg] [options]'
            result << ''
            result << 'Commands:'
            commands.collect{ it.name }.sort().each { result << "  $it".toString()  }
            result << ''
        }
        else {
            def sub = commands.find { it.name == args[0] }
            if( sub )
                sub.usage(result)
            else {
                throw new AbortOperationException("Unknown cloud sub-command: ${args[0]}")
            }
        }

        println result.join('\n').toString()
    }


    private void addField(String fieldName, List<String> result) {
        def annot = this.class.getDeclaredField(fieldName)?.getAnnotation(Parameter)
        if( annot ) {
            result << '  ' + annot.names().join(', ')
            result << '     ' + annot.description()
        }
        else {
            log.debug "Unknown help field: $fieldName"
        }
    }

    private void launchCommonOpts(List<String> result) {
        result << 'Options:'
        addField('instanceCount', result)
        addField('driver', result)
        addField('imageId', result)
        addField('instanceType', result)
        addField('spotPrice', result)
        addField('region', result)
        addField('profile', result)
    }

    /**
     * Implements the cloud `list-node` sub-command
     */
    class ListNodes implements SubCmd {

        @Override
        String getName() { 'list' }

        String getShortcut() { 'l' }

        @Override
        void apply(String clusterName) {

            if( !clusterName ) {
                printAvailableClusterNames()
            }
            else {
                printClusterMembers(clusterName)
            }

        }

        private printClusterMembers(String clusterName) {
            def tags = [:]; tags[TAG_CLUSTER_NAME] = clusterName
            def builder = new TableBuilder()
                .head('INSTANCE ID')
                .head('ADDRESS')
                .head('STATUS')
                .head('ROLE')

            driver.eachInstanceWithTags(tags) { CloudInstance item ->
                builder << item.id
                builder << item.address
                builder << item.state
                builder << item.role
                builder.closeRow()
            }

            println builder.toString()
        }

        private printAvailableClusterNames() {
            def names = getClusterNames()

            if( !names ) {
                println "No cluster available"
            }
            else {
                println "The following clusters are available: "
                names.sort().each { println "  $it"  }
            }
        }

        @Override
        void usage(List<String> result) {
            result << 'List cloud cluster node(s)'
            result << "Usage: nextflow cloud $name <cluster name> [options]".toString()
            result << ''
            addField('region', result)
            result << ''
        }
    }

    /**
     * Implements the cloud `launch-cluster` sub-command
     */
    class CreateCluster implements SubCmd {

        @Override
        String getName() { 'create' }

        String getShortcut() { 'c' }

        @Override
        void apply(String clusterName) {
            checkName(clusterName)
            int count = instanceCount ?: 1

            // -- show the config and ask for confirmation
            def cloudConfig = makeConfig(clusterName)
            printConfig(cloudConfig, count)
            if( !yes )
                prompt("Please confirm you really want to launch the cluster with above configuration")

            // -- make sure the name it's not already used
            if( clusterName in getClusterNames() )
                throw new AbortOperationException("A cluster with name `$clusterName` already exists -- Cannot continue")

            if( count > 1 ) {
                // -- launch the worker nodes
                def workerIds = addNodes0((int)count-1, cloudConfig.setRole(ROLE_WORKER))
                printWorkerInstances(workerIds)
            }

            // -- launch the master node
            def masterId = addNodes0(1, cloudConfig.setRole(ROLE_MASTER))
            printMasterInstance(masterId.first(), cloudConfig)
        }

        @Override
        void usage(List<String> result) {
            result << 'Deploy a Nextflow cluster in the cloud'
            result << "Usage: nextflow cloud $name <cluster name> [options]".toString()
            result << ''
            launchCommonOpts(result)
            addField('yes', result)
            result << ''
        }
    }


    /**
     * Implements the cloud `launch-nodes` sub-command
     */
    class JoinCluster implements SubCmd {

        @Override
        String getName() { return 'join' }

        String getShortcut() { 'j' }

        @Override
        void apply(String clusterName) {
            checkName(clusterName)
            int count = instanceCount ?: 1

            // -- show the config and ask for confirmation
            def cloudConfig = makeConfig(clusterName).setRole(ROLE_WORKER)
            printConfig(cloudConfig, count)
            prompt("Please confirm you really want to launch the cluster nodes(s) with the above configuration: `$clusterName`")

            if( !getClusterNames().contains(clusterName) )
                throw new AbortOperationException("No cluster available with name `$clusterName`")

            // launch the node
            def instanceIds = addNodes0( count, cloudConfig)
            printWorkerInstances(instanceIds)
        }

        @Override
        void usage(List<String> result) {
            result << 'Add one or more nodes to a running cluster'
            result << "Usage: nextflow cloud $name <cluster name> [options]".toString()
            result << ''
            launchCommonOpts(result)
            result << ''
        }
    }

    /**
     * Implements the cloud `spot-prices` sub-command
     */
    class SpotPrices implements SubCmd {

        @Override
        String getName() { return 'spot-prices' }

        String getShortcut() { 'sp' }


        @Override
        void apply(String none) {
            final all = CmdCloud.this.all
            final history = CmdCloud.this.history
            final filter = history ? "type=='$history'" : fixEqualsOp(CmdCloud.this.filter)
            final sort = CmdCloud.this.sort ?: 'price'
            final Script filterScript = filter ? new Grengine().create("{ -> $filter }") : null
            log.debug "Spot-price: filter=`$profile`; sort=`$sort`; all=$all; history=`$history`"

            def allPrices = new ArrayList<Map>()
            driver.eachSpotPrice { CloudSpotPrice entry ->

                if( entry.description.contains('Windows') ) {
                    return  // ignore Windows entries
                }

                def type = driver.describeInstanceType(entry.type)
                if( !type ) {
                    log.warn1 "Unknown instance type: ${entry.type}"
                    return
                }

                def record = [
                        type: entry.type,
                        zone: entry.zone,
                        price: entry.price as float,
                        pricecpu: Math.round((entry.price as float) * 10_000 / type.cpus) / 10_000,
                        description: entry.description,
                        cpus:  type.cpus,
                        mem: type.memory,
                        disk: type.disk,
                        numOfDisks: type.numOfDisks,
                        time: entry.timestamp
                ]

                def accept = true
                if( filterScript ) {
                    def binding = new Binding(record)
                    filterScript.setBinding(binding)
                    accept = ((Closure)filterScript.run()).call()
                }

                if( accept ) {
                    allPrices << record
                }
            }
            log.debug "Found spot-price records: ${allPrices.size()}"


            if( all ) {
                printAll(allPrices, sort)
            }
            else if( history ) {
                printHistory(allPrices)
            }
            else {
                printCheapest(allPrices, sort)
            }
        }

        private Comparator<Map> compareBy( String sort, List<Map> records=null ) {

            if( !sort ) return null

            // -- multi-column comparator
            def fields = sort.tokenize(', ');

            // validate fields
            if( records ) {
                def allFields = records[0].keySet()
                def copy = new ArrayList<>(fields)
                copy.removeAll(allFields)
                if( copy ) {
                    throw new AbortOperationException("Not a valid sort field(s): ${copy.join(',')}")
                }
            }

            return { Map o1, Map o2 ->
                for( int i=0; i<fields.size(); i++ )  {
                    def f = fields[i]
                    def v = o1.get(f) <=> o2.get(f)
                    if( v != 0 ) return v
                }
                return 0
            } as Comparator<Map>

        }

        private void printAll(List<Map> allPrices, String sort) {

            def prices = sort ? allPrices.sort(false,compareBy(sort, allPrices)) : allPrices

            def table = new TableBuilder()
                    .head('TYPE')
                    .head('PRICE')
                    .head('PRICE/CPU')
                    .head('ZONE')
                    .head('DESCRIPTION')
                    .head('CPUS')
                    .head('MEMORY', TextLabel.Align.RIGHT)
                    .head('DISK')
                    .head('TIMESTAMP')

            // show
            prices.each { entry ->
                table << entry.type
                table << String.format("%.4f", entry.price)
                table << String.format("%.4f", entry.pricecpu)
                table << entry.zone
                table << entry.description
                table << entry.cpus
                table << entry.mem
                table << (entry.disk ? "${entry.numOfDisks} x ${entry.disk}" : '-')
                table << SysHelper.fmtDate((Date)entry.time)
                table.closeRow()
            }

            println table.toString()
        }

        private void printCheapest(List<Map> allPrices, String sort) {

            def filter = new HashMap<String,Map>()
            def itr = allPrices.sort(false,compareBy('type,time', allPrices)).iterator()
            while( itr.hasNext() ) {
                def rec = itr.next()
                //println rec

                final id = (String)rec.type
                final prev = filter[id]
                if( !prev ) {
                    filter[id] = rec
                    continue
                }

                if( rec.zone == prev.zone) {
                    // it's a most recent rec update it, take it
                    filter[id] = rec
                }
                else if( rec.price < prev.price ) {
                    // it's a cheapest price in a different zone, take it
                    filter[id] = rec
                }
            }

            // order with the sorting criteria provided by the user
            def prices = sort ? filter.values().sort(false,compareBy(sort,allPrices)) : filter.values()

            // render the table
            def table = new TableBuilder()
                    .head('TYPE')
                    .head('PRICE')
                    .head('PRICE/CPU')
                    .head('ZONE')
                    .head('DESCRIPTION')
                    .head('CPUS')
                    .head('MEMORY', TextLabel.Align.RIGHT)
                    .head('DISK')

            // show
            prices.each { entry ->
                table << entry.type
                table << String.format("%.4f", entry.price)
                table << String.format("%.4f", entry.pricecpu)
                table << entry.zone
                table << entry.description
                table << entry.cpus
                table << entry.mem
                table << (entry.disk ? "${entry.numOfDisks} x ${entry.disk}" : '-')
                table.closeRow()
            }

            println table.toString()
        }

        private void printHistory(List<Map> allPrices) {

            allPrices.sort(true, compareBy('zone,time',allPrices))

            def table = new TableBuilder()
                    .head('ZONE')
                    .head('TIMESTAMP')
                    .head('PRICE')
                    .head('PRICE/CPU')

            def z = null
            allPrices.each { entry ->
                if( z != entry.zone ) {
                    z = entry.zone
                    table << z
                    table.closeRow()
                }

                table << ''
                table << SysHelper.fmtDate((Date)entry.time)
                table << String.format("%.4f", entry.price)
                table << String.format("%.4f", entry.pricecpu)
                table.closeRow()
            }

            println table.toString()
        }

        @Override
        void usage(List<String> result) {
            result << 'Show current spot instance prices'
            result << "Usage: nextflow cloud $name [filter]".toString()
            result << ''
            result << 'Options:'
            addField('all', result)
            addField('filter', result)
            addField('history', result)
            addField('sort', result)
            addField('region', result)
            result << ''
        }
    }

    /**
     * Implements the cloud `shutdown` sub-command
     */
    class Shutdown implements SubCmd {

        @Override
        String getName() { return 'shutdown' }

        String getShortcut() { null }

        @Override
        void apply(String clusterName) {
            checkName(clusterName)

            def instanceIds = []
            def tags = [:]; tags[TAG_CLUSTER_NAME] = clusterName

            driver.eachInstanceWithTags(tags) { CloudInstance it -> instanceIds.add(it.id) }
            if( !instanceIds ) {
                println "No cluster found for name: `$clusterName`"
                return
            }

            prompt("Please confirm you really want to shutdown cluster: `$clusterName`")
            driver.terminateInstances(instanceIds)
        }

        @Override
        void usage(List<String> result) {
            result << 'Shutdown a running cluster the cloud'
            result << "Usage: nextflow cloud $name <cluster name> [options]".toString()
            result << ''
            result << 'Options:'
            addField('region', result)
            result << ''
        }
    }
}
