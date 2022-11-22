/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.script

import groovy.transform.CompileStatic

/**
 * Placeholder for process inputs
 *
 * @author Jacob E Munro <jacob.e.munro@gmail.com>
 */

@CompileStatic
class InputPlaceholder {

    private Integer index

    private String name

    InputPlaceholder() { }

    InputPlaceholder (Integer index) {
        this.index = index
    }

    InputPlaceholder (String name) {
        this.name = name
    }

    InputPlaceholder getAt(Integer i) {
        new InputPlaceholder(i)
    }

    InputPlaceholder getAt(String name) {
        new InputPlaceholder(name)
    }

    @Override
    Object getProperty(String name) {
        this.hasProperty(name) ? this.getProperties()[name] : new InputPlaceholder(name)
    }

    Object getInput(Object object) {
        
        if (object instanceof ChannelOut) {
            ChannelOut channelOut = object as ChannelOut
            return name ? channelOut.getProperty(name) : index ? channelOut[index] : channelOut[0]
        }
        return object
    }

    String getName() { name }
}