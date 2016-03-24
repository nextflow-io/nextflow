/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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


import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter;

public class JCommanderExample {

    @Parameter
    private List<String> parameters = new ArrayList<String>();

    @Parameter(names = [ "-log", "-verbose" ], description = "Level of verbosity")
    private Integer verbose = 1;

    @Parameter(names = "-groups", description = "Comma-separated list of group names to be run")
    private String groups;

    @Parameter(names = "-debug", description = "Debug mode")
    private boolean debug = false;

    @DynamicParameter(names = "--\$", description = "Dynamic parameters go here", assignment = ' ')
    private Map<String, String> params = new HashMap<String, String>();

}



JCommanderExample jct = new JCommanderExample();
String[] argv = ["-log", "2", "-groups", "unit", 'a', '--$b','c']
new JCommander(jct, argv);

assert jct.verbose == 2
assert jct.parameters == ['a']
assert jct.params == [b:'c']