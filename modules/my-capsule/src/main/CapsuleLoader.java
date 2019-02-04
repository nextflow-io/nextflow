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

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Implement capsule custom strategy to add maven dependencies as
 * a command line argument
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CapsuleLoader extends Capsule {

    /**
     * Constructs a capsule from the given JAR file
     *
     * @param jarFile  the path to the JAR file
     */
    protected CapsuleLoader(Path jarFile) {
        super(jarFile);
    }

    @Override
    protected ProcessBuilder prelaunch(List<String> jvmArgs, List<String> args) {
        ProcessBuilder pb = super.prelaunch(jvmArgs, args);

        String drip = System.getenv().get("NXF_DRIP");
        if( drip != null && !drip.isEmpty() ) {
            pb.command().set(0, drip);
            return pb;
        }

        return pb;
    }

    @Override
    @SuppressWarnings("unchecked")
    protected <T> T attribute(Map.Entry<String, T> attr) {

        if (ATTR_DEPENDENCIES == attr && System.getenv("NXF_GRAB") != null ) {
            String deps = System.getenv("NXF_GRAB");
            List<String> parent = super.attribute((Map.Entry<String,List<String>>)attr);
            return (T)extendDepsWith(deps, parent);
        }

        if (ATTR_APP_CLASS_PATH == attr && System.getenv("NXF_CLASSPATH") != null) {
            String classpath = System.getenv("NXF_CLASSPATH");
            List<String> parent = super.attribute((Map.Entry<String, List<String>>)attr);
            return (T)extendClassPathWith(classpath, parent);
        }

        return super.attribute(attr);
    }

    static List<String> extendClassPathWith(String classpath, List<String> origin) {
        if( classpath == null || "".equals(classpath.trim()) ) {
            return origin;
        }

        List<String> result = new ArrayList<>();
        if( origin != null ) {
            result.addAll(origin);
        }

        for( String lib : classpath.split(":") ) {
            String trimmed = lib.trim();
            if( trimmed.length()>0 )
                result.add(trimmed);
        }

        return result;
    }

    /*
     * get a blank separated list of dependencies
     * having the "maven" syntax i.e. group:name:version
     */
    static List<String> extendDepsWith(String dependencies, List<String> origin) {
        if( dependencies == null || "".equals(dependencies.trim()) ) {
            return origin;
        }

        // split by black
        List<String> result = new ArrayList<>();
        for( String dep : dependencies.split(" |\n") ) {
            String trimmed = dep.trim();
            if( trimmed.length()>0 )
                result.add(trimmed);
        }

        // add all the inherited if any
        if( origin != null && origin.size()>0 ) {
            result.addAll(0, origin);
        }

        return result;
    }

}
