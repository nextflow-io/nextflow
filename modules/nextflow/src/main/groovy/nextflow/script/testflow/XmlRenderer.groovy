package nextflow.script.testflow

import javax.xml.stream.XMLOutputFactory
import javax.xml.stream.XMLStreamWriter
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.time.Duration

/**
 *
 * To create a test report XML file using XUnit format (compatible to JUnit format)
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
class XmlRenderer {

    /**
     * Write XML testsuite file
     *
     * @param suite Testsuite to write
     * @param testDir Output folder
     */
    static void write(TestSuite suite, Path testDir) {
        final xmlFile = testDir.resolve("TEST-${suite.name}.xml")
        final Charset charset = Charset.defaultCharset()
        testDir.mkdirs()
        Writer bw = Files.newBufferedWriter(xmlFile, charset)
        final XMLOutputFactory xof = XMLOutputFactory.newFactory()
        final XMLStreamWriter w = xof.createXMLStreamWriter(bw)
        w.writeStartDocument(charset.displayName(), "1.0")
        renderTestSuite(w, suite)
        w.writeEndDocument()
        w.flush()
        bw.flush()
        bw.close()
    }

    static private void renderTestSuite(XMLStreamWriter w, TestSuite s) {
        // testsuite
        w.writeStartElement("testsuite")
        w.writeAttribute("name", s.name)
        w.writeAttribute("tests", String.format("%d", s.tests))
        w.writeAttribute("skipped", String.format("%d", s.skipped))
        w.writeAttribute("failures", String.format("%d", s.failures))
        w.writeAttribute("errors", String.format("%d", s.errors))
        w.writeAttribute("timestamp", s.timestamp.toString())
        w.writeAttribute("time", formatTime(s.time))

        // testcase
        s.testcase.each { renderTestCase(w, it) }

        // system-out
        w.writeStartElement("system-out")
        w.writeCData(s.systemOut)
        w.writeEndElement()

        // system-err
        w.writeStartElement("system-err")
        w.writeCData(s.systemErr)
        w.writeEndElement()

        w.writeEndElement()
    }

    static private void renderTestCase(XMLStreamWriter w, TestCase t) {

        w.writeStartElement("testcase")
        w.writeAttribute("name", t.name)
        w.writeAttribute("classname", t.className)
        w.writeAttribute("time", formatTime(t.time))

        if (t.failure) {
            w.writeStartElement("failure")
            w.writeAttribute("message", t.failure.message)
            w.writeAttribute("type", t.failure.type)
            w.writeCData(t.failure.content)
            w.writeEndElement()
        }

        w.writeEndElement()
    }

    static private String formatTime(Duration d) {
        String.format("%.4f", d.toMillis() / 1000.0f)
    }
}
