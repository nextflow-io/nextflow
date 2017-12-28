/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
