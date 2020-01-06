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

package nextflow.dag
import java.nio.charset.Charset
import java.nio.file.Path
import java.text.SimpleDateFormat;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

/**
 * Render the DAG in Gexf/Gephi
 * See https://gephi.org/ for more info.
 *
 * @author Pierre Lindenbaum / yokofakun
 */
class GexfRenderer implements DagRenderer {
    private static final String XMLNS="http://www.gexf.net/1.3";
    private static final String XSI_SCHEMA_LOCATION="http://www.gexf.net/1.3 http://www.gexf.net/1.3/gexf.xsd";
    private static final String VERSION="1.3";


    private final String name;
    GexfRenderer(String name) {
    this.name = name;    
    }


    @Override
    void renderDocument(DAG dag, Path file) {
    final Charset charset = Charset.defaultCharset();
    Writer bw = Files.newBufferedWriter(file, charset); 
    final XMLOutputFactory xof = XMLOutputFactory.newFactory();
    final XMLStreamWriter w = xof.createXMLStreamWriter(bw);    
    w.writeStartDocument(charset.displayName(),"1.0");
    w.writeStartElement("gexf");
    w.writeAttribute("xmlns",XMLNS);
    w.writeAttribute("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance");
    w.writeAttribute("xsi:schemaLocation",XSI_SCHEMA_LOCATION);
    w.writeAttribute("version", VERSION);

    /* meta */
    w.writeStartElement("meta");
        w.writeAttribute("lastmodifieddate",new SimpleDateFormat("yyyy-MM-dd").format(new Date()));
        w.writeStartElement("creator");
          w.writeCharacters("Nextflow");
        w.writeEndElement();
        w.writeStartElement("description");
          w.writeCharacters(this.name);
        w.writeEndElement();
    w.writeEndElement();
            
    /* graph */
    w.writeStartElement("graph");
    w.writeAttribute("mode", "static");
    w.writeAttribute("defaultedgetype", "directed");
    
    
    /* attributes */
    w.writeStartElement("attributes");
    w.writeAttribute("class","node");
    w.writeAttribute("mode","static");
    w.writeEmptyElement("attribute");
        w.writeAttribute("id", "type");
        w.writeAttribute("title","type");
        w.writeAttribute("type", "string");
    w.writeEndElement();//attributes
    
            
    /* vertex/node */
    w.writeStartElement("nodes");
    dag.vertices.each { vertex -> renderVertex(w, vertex ) }
    w.writeEndElement();

    /* edges */
    w.writeStartElement("edges");
    dag.edges.each { edge -> renderEdge(w, edge ) }
    w.writeEndElement();
            
    w.writeEndElement();
    w.writeEndDocument();
    w.flush();
    bw.flush();
    bw.close();
    }

    private void renderVertex(w,vertex) {

    w.writeStartElement("node");
    w.writeAttribute("id",vertex.getName());
    w.writeAttribute("label",vertex.label?vertex.label:vertex.getName());
    
    w.writeStartElement("attvalues");
        w.writeEmptyElement("attvalue");
            w.writeAttribute("for","type");
            w.writeAttribute("value",vertex.type);
    w.writeEndElement();//attvalues

    w.writeEndElement();//node
    }


    private void renderEdge(w,edge) {
        assert edge.from != null && edge.to != null
    w.writeStartElement("edge");
    w.writeAttribute("type", "directed");
    w.writeAttribute("source",edge.from.name);
    w.writeAttribute("target",edge.to.name);
    if(edge.label) w.writeAttribute("label",edge.label);
    w.writeEndElement();//edge
    }

}
