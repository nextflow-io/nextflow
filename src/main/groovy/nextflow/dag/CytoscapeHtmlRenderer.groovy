package nextflow.dag
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

import groovy.transform.CompileStatic
import java.nio.file.Path

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class CytoscapeHtmlRenderer implements DagRenderer {

    /**
     * Render the DAG in HTML using Cytoscape.js
     * to the specified file.
     * See http://js.cytoscape.org for more info.
     */
    @Override
    void renderDocument(DAG dag, Path file) {
        String tmplPage = readTemplate()
        String network = CytoscapeJsRenderer.renderNetwork(dag)
        file.text = tmplPage.replaceAll(~/\/\* REPLACE_WITH_NETWORK_DATA \*\//, network)
    }

    private String readTemplate() {
        StringWriter writer = new StringWriter();
        def res = CytoscapeJsRenderer.class.getResourceAsStream('cytoscape.js.dag.template.html')
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }
}
