/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.plugin.hello

import java.util.concurrent.CompletableFuture

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Global
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import io.nextflow.gradle.extensions.Function
import io.nextflow.gradle.extensions.Operator
import io.nextflow.gradle.extensions.PluginExtensionPoint

/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Slf4j
class HelloExtension extends PluginExtensionPoint {

    /*
    * Flag to check if init was called
     */
    boolean initialized = false

    /*
     * A session hold information about current execution of the script
     */
    private Session session

    /*
     * nf-core initializes the plugin once loaded and session is ready
     * @param session
     */
    @Override
    protected void init(Object session) {
        this.session = session as Session
        this.initialized = true
    }

    /*
     * reverse is a `producer` method and will be available to the script because:
     *
     * - it's public
     * - it returns a DataflowWriteChannel
     *
     * nf-core will inspect the extension class and allow the script to call all these kind of methods
     *
     * the method can require arguments but it's not mandatory, it depends of the business logic of the method
     *
     * business logic can write into the channel once ready and values will be consumed from it
     */
    DataflowWriteChannel reverse(String message) {
        createReverseChannel(message)
    }

    /*
     * this Factory can't be imported because has not the right signature
     */
    String reverseCantBeImportedBecauseWrongSignature(String message){
        message
    }

    /*
    * goodbye is a `consumer` method as it receives values from a channel to perform some logic.
    *
    * Consumer methods are introspected by nextflow-core and include into the DSL if the method:
    *
    * - it's public
    * - it returns a DataflowWriteChannel
    * - it has only one arguments of DataflowReadChannel class
    *
    * a consumer method needs to proporcionate 2 closures:
    * - a closure to consume items (one by one)
    * - a finalizer closure
    *
    * in this case `goodbye` will consume a message and will store it as an upper case
    */
    @Operator
    DataflowWriteChannel goodbye(DataflowReadChannel source) {
        final target = CH.createBy(source)
        final next = {
            target.bind(it)
        }
        final done = {
            target.bind(Channel.STOP)
        }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    }

    /*
    * this Operator can't be imported because has not the right signature
     */
    @Operator
    String goodbyeWrongSignature(DataflowReadChannel source) {
        source.val
    }

    protected DataflowWriteChannel createReverseChannel(final String message){
        final channel = CH.create()
        if( NF.isDsl2() ){
            session.addIgniter { ->
                businessLogicHere(channel, message)
            }
        }else{
            businessLogicHere(channel, message)
        }
        channel
    }

    HelloFunctions functions = new HelloFunctions()

    /**
     * An annotate @Function function to be imported as a custom function by the parser
     * @param lang
     * @return
     */
    @Function
    String sayHello(String lang='en'){
        assert initialized, "PluginExtension was not initialized"
        // sayHello is the entrypoint where we can write all the logic or delegate to other classes, ...
        return functions.sayHello(lang)
    }

    String aNonImportedFunction(){
        throw new IllegalAccessException("This function can't be imported")
    }

    /*
    * businessLogicHere will send, across the channel, the message reversed
    * and after will send an STOP signal to let know the channel it has been finished
    */
    protected static businessLogicHere(final DataflowWriteChannel channel, final String message){
        def future = CompletableFuture.runAsync({
            channel.bind(message.reverse())
            channel.bind(Channel.STOP)
        })
        future.exceptionally(this.&handlerException)
    }

    /*
    * an util class to trace exceptions
    */
    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }
}
