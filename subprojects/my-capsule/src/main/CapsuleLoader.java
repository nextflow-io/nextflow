/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

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
     * @param cacheDir the path to the (shared) Capsule cache directory
     */
    protected CapsuleLoader(Path jarFile, Path cacheDir) {
        super(jarFile, cacheDir);
    }


    protected ProcessBuilder prelaunch(List<String> jvmArgs, List<String> args) {
        ProcessBuilder pb = super.prelaunch(jvmArgs,args);

        String drip = System.getenv().get("NXF_DRIP");
        if( drip != null && !"".equals(drip) ) {
            pb.command().set(0, drip);
            return pb;
        }

        String javaCmd = System.getenv().get("JAVA_CMD");
        if( javaCmd != null || !"".equals(javaCmd) ) {
            // use the Java command provided by the env variable
            pb.command().set(0, javaCmd);
        }

        return pb;
    }

    @Override
    protected List<String> getDependencies() {
        List<String> parent = super.getDependencies();
        return extendDepsWith(System.getenv("NXF_GRAB"), parent);
    }

    protected List<Path> buildClassPath() {
        List<Path> parent = super.buildClassPath();
        String classpath = System.getenv("NXF_CLASSPATH");
        return extendClassPathWith(classpath, parent);
    }

    static List<Path> extendClassPathWith(String classpath, List<Path> origin) {
        if( classpath == null || "".equals(classpath.trim()) ) {
            return origin;
        }

        List<Path> result = new ArrayList<>();
        if( origin != null ) {
            result.addAll(origin);
        }


        for( String lib : classpath.split(":") ) {
            String trimmed = lib.trim();
            if( trimmed.length()>0 )
                result.add(Paths.get(trimmed));
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
