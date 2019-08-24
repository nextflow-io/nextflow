/*
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


/**
 * Holds process script meta-data
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BodyDef implements Cloneable {

    /** The actual task code */
    Closure closure

    /** The task code as source, eventually used to report diagnostic message in case of error */
    String source

    /** The script type, either a system text script or groovy native code*/
    ScriptType type

    /** When true the script is interpreted a shell script in which $xxx variable are not interpreted by groovy */
    boolean isShell

    /**
     * The set of token variable and properties referenced in the task script
     */
    Set<TokenValRef> valRefs

    /**
     * Hold the use provided task body
     *
     * @param closure The task body closure
     * @param source The task source string
     * @param section The section in which it has been, declared. One of
     *      {@code exec}, {@code script} or {@code shell} literals
     */
    BodyDef(Closure closure, String source, String section = 'script' ) {
        this.closure = closure
        this.source = source
        this.valRefs = Collections.emptySet()
        setType(section)
    }

    BodyDef(Closure closure, String source, String section, TokenValRef... values ) {
        this(closure, source, section)
        this.valRefs = values != null ? values as Set : Collections.<TokenValRef>emptySet()
    }

    BodyDef clone() {
        def result = (BodyDef) super.clone()
        result.closure = (Closure) closure.clone()
        result.valRefs = new HashSet(valRefs)
        return result
    }
    /**
     * Given a code section literal sets the task body {@link #type} and {@link #isShell} flag
     *
     * @param section
     * @return
     */
    protected BodyDef setType(String section )  {

        switch( section ) {
            case 'exec':
                type = ScriptType.GROOVY
                isShell = false
                break

            case 'script':
                type = ScriptType.SCRIPTLET
                isShell = false
                break

            case 'shell':
                this.type = ScriptType.SCRIPTLET
                isShell = true
                break

            case 'workflow':
                type = ScriptType.GROOVY
                isShell = false
                break

            default:
                throw new IllegalArgumentException("Not a valid process body definition: $section")
        }

        return this
    }

    String toString() {
        "BodyDef[closure: $closure; types: $type; shell: $isShell; variables: $valRefs; source:\n$source \n]"
    }

    /**
     * @return
     *      The set of variable and properties referenced in the user script.
     *      NOTE: it includes properties in the form {@code object.propertyName}
     */
    Set<String> getValNames() {
        valRefs *. name
    }
}
