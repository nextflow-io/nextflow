package nextflow.script.testflow

import groovy.xml.MarkupBuilder
import nextflow.Const

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.time.Duration
import java.time.LocalDateTime

/**
 * Creates an HTML test report for a collection of test suites
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
class HtmlRenderer {

    /**
     * Write test HTML report
     *
     * @param testRun A test run of testSuites to report
     * @param reportDir Output folder to write the HTML
     */
    static void write(TestRun testRun, Path reportDir) {
        reportDir.mkdirs()
        writeCss(reportDir)
        writeJs(reportDir)
        testRun.suites.each { writeTestSuite(it, reportDir) }
        writeIndex(testRun, reportDir)
    }

    private static void writeIndex(TestRun testRun, Path reportDir) {
        final indexFile = reportDir.resolve("index.html")
        final Charset charset = Charset.defaultCharset()
        Writer writer = Files.newBufferedWriter(indexFile, charset)

        writer.write("<!DOCTYPE html>\n")
        def builder = new MarkupBuilder(writer)
        builder.html {
            head {
                meta('http-equiv': 'Content-Type', content: 'text/html; charset=utf-8')
                meta('http-equiv': 'x-ua-compatible', content: 'IE=edge')
                title("Test results - Test Summary")
                link(href: 'style.css', rel: 'stylesheet', type: 'text/css')
                script('', type: 'text/javascript', src: 'report.js')
            }

            body {
                div(id: 'content') {
                    h1("Test Summary")
                    div(id: 'summary') {
                        table {
                            tr {
                                td {
                                    div(class: 'summaryGroup') {
                                        table {
                                            tr {
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(testRun.completed, class: 'counter')
                                                        p("tests")
                                                    }
                                                }
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(testRun.failed, class: 'counter')
                                                        p("failures")
                                                    }
                                                }
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(testRun.skipped, class: 'counter')
                                                        p("ignored")
                                                    }
                                                }
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(humanReadable(testRun.duration), class: 'counter')
                                                        p("duration")
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                td {
                                    div(class: 'infoBox' + (testRun.failed ? ' failures' : ''), id: 'successRate') {
                                        div("${testRun.successRate}%", class: 'percent')
                                        p("successful")
                                    }
                                }
                            }
                        }
                    }

                    div(id: 'tabs') {
                        ul(class: 'tabLinks') {
                            li { a("Failed tests", href: '#tab0') }
                            li { a("Tests", href: '#tab1') }
                        }

                        div(id: 'tab0', class: 'tab') {
                            h2("Failed tests")
                            ul(class: 'linkList') {
                                testRun.suites.each { suite ->
                                    suite.testcase.each { testcase ->
                                        testcase.failed ? li {
                                            a(suite.name, href: "${suite.name}.html")
                                            span(" > ")
                                            a(testcase.name, href: "${suite.name}.html#${testcase.name}")
                                        } : null
                                    }
                                }
                            }
                        }

                        div(id: 'tab1', class: 'tab') {
                            h2("Tests")
                            table {
                                thead {
                                    tr {
                                        th("Script")
                                        th("Tests")
                                        th("Failures")
                                        th("Ignored")
                                        th("Duration")
                                        th("success rate")
                                    }
                                }
                                tbody {
                                    testRun.suites.each {testsuite ->
                                       tr {
                                           td(class:"failures") {
                                               a(testsuite.name, href:"${testsuite.name}.html")
                                           }
                                           td(testsuite.tests)
                                           td(testsuite.failures + testsuite.errors)
                                           td(testsuite.skipped)
                                           td(humanReadable(testsuite.time))
                                           td("${((testsuite.tests - testsuite.errors - testsuite.failures) * 100 / testsuite.tests)}%", class:'failures')
                                       }
                                    }
                                }
                            }
                        }
                    }

                    div(id: 'footer') {
                        p("Generated by Nextflow ${Const.APP_VER} at ${LocalDateTime.now()}")
                    }
                }
            }
        }

        writer.flush()
        writer.close()
    }

    private static void writeTestSuite(TestSuite testSuite, Path reportDir) {

        final completed = testSuite.tests
        final skipped = testSuite.skipped
        final failures = testSuite.failures
        final errors = testSuite.errors
        final duration = testSuite.time
        final failed = failures + errors
        final successRate = ((completed - failed) * 100 / completed) as int

        final outFile = reportDir.resolve("${testSuite.name}.html")
        final Charset charset = Charset.defaultCharset()
        Writer writer = Files.newBufferedWriter(outFile, charset)
        writer.write("<!DOCTYPE html>\n")
        def builder = new MarkupBuilder(writer)
        builder.omitNullAttributes = true
        builder.html {
            head {
                meta('http-equiv': 'Content-Type', content: 'text/html; charset=utf-8')
                meta('http-equiv': 'x-ua-compatible', content: 'IE=edge')
                title("Test results - ${testSuite.name}")
                link(href: 'style.css', rel: 'stylesheet', type: 'text/css')
                script('', type: 'text/javascript', src: 'report.js')
            }
            body {
                div(id: 'content') {
                    h1("Test ${testSuite.name}")
                    div(class: 'breadcrumbs') {
                        a("all", href: 'index.html')
                        span(" > ${testSuite.name}")
                        a("[workdir]",  href: relativeURL(reportDir, testSuite.workDir))
                    }
                    div(id: 'summary') {
                        table {
                            tr {
                                td {
                                    div(class: 'summaryGroup') {
                                        table {
                                            tr {
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(completed, class: 'counter')
                                                        p("tests")
                                                    }
                                                }
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(failed, class: 'counter')
                                                        p("failures")
                                                    }
                                                }
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(skipped, class: 'counter')
                                                        p("ignored")
                                                    }
                                                }
                                                td {
                                                    div(class: 'infoBox', id: 'tests') {
                                                        div(humanReadable(duration), class: 'counter')
                                                        p("duration")
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                td {
                                    div(class: 'infoBox' + (failed > 0 ? ' failures' : ''), id: 'successRate') {
                                        div("${successRate}%", class: 'percent')
                                        p("successful")
                                    }
                                }
                            }
                        }
                    }

                    div(id: 'tabs') {
                        ul(class: 'tabLinks') {
                            li { a("Failed tests", href: '#tab0') }
                            li { a("Tests", href: '#tab1') }
                        }

                        div(id: 'tab0', class: 'tab') {
                            h2("Failed tests")
                            testSuite.testcase.each { testcase ->
                                testcase.failed ? div(class: 'test') {
                                    a("", name: testcase.name)
                                    h3(testcase.name, class: 'failures')
                                    testcase.workDir ? a("[${processDirReference(testcase.workDir)}]", href:relativeURL(reportDir, testcase.workDir), class:'workdir') : null
                                    span(class: 'code') { pre(testcase.failure.content) }
                                } : null
                            }
                        }

                        div(id: 'tab1', class: 'tab') {
                            h2("Tests")
                            table {
                                thead {
                                    tr {
                                        th("Test")
                                        th("Duration")
                                        th("Result")
                                    }
                                }
                                tbody {
                                    testSuite.testcase.each { testcase ->
                                        tr {
                                            td(class: testcase.failed ? "failures" : "success") {
                                                span(testcase.name)
                                                testcase.workDir ? a("[${processDirReference(testcase.workDir)}]", href:relativeURL(reportDir, testcase.workDir), class:'workdir') : null
                                            }
                                            td(humanReadable(testcase.time), class: testcase.failed ? "failures" : "success")
                                            td(testcase.failed ? "failed" : "success", class: testcase.failed ? "failures" : "success")
                                        }
                                    }
                                }
                            }
                        }
                    }

                    div(id: 'footer') {
                        p("Generated by Nextflow ${Const.APP_VER} at ${LocalDateTime.now()}")
                    }
                }
            }
        }

        writer.flush()
        writer.close()
    }

    static String humanReadable(Duration duration) {
        return duration.toString()
                .substring(2)
                .replaceAll('(\\d[HMS])(?!$)', '$1 ')
                .toLowerCase()
    }

    static String relativeURL(Path base, Path inner) {
        base.toAbsolutePath().relativize(inner)
    }

    static String processDirReference(Path processDir) {
        "${processDir.parent.name}/${processDir.name.substring(0, 6)}"
    }

    private static void writeJs(Path reportDir) {
        final indexFile = reportDir.resolve("report.js")
        final Charset charset = Charset.defaultCharset()
        Writer writer = Files.newBufferedWriter(indexFile, charset)

        writer.write("""
(function (window, document) {
    "use strict";

    var tabs = {};

    function changeElementClass(element, classValue) {
        if (element.getAttribute("className")) {
            element.setAttribute("className", classValue);
        } else {
            element.setAttribute("class", classValue);
        }
    }

    function getClassAttribute(element) {
        if (element.getAttribute("className")) {
            return element.getAttribute("className");
        } else {
            return element.getAttribute("class");
        }
    }

    function addClass(element, classValue) {
        changeElementClass(element, getClassAttribute(element) + " " + classValue);
    }

    function removeClass(element, classValue) {
        changeElementClass(element, getClassAttribute(element).replace(classValue, ""));
    }

    function initTabs() {
        var container = document.getElementById("tabs");

        tabs.tabs = findTabs(container);
        tabs.titles = findTitles(tabs.tabs);
        tabs.headers = findHeaders(container);
        tabs.select = select;
        tabs.deselectAll = deselectAll;
        tabs.select(0);

        return true;
    }

    function getCheckBox() {
        return document.getElementById("line-wrapping-toggle");
    }

    function getLabelForCheckBox() {
        return document.getElementById("label-for-line-wrapping-toggle");
    }

    function findCodeBlocks() {
        var spans = document.getElementById("tabs").getElementsByTagName("span");
        var codeBlocks = [];
        for (var i = 0; i < spans.length; ++i) {
            if (spans[i].className.indexOf("code") >= 0) {
                codeBlocks.push(spans[i]);
            }
        }
        return codeBlocks;
    }

    function forAllCodeBlocks(operation) {
        var codeBlocks = findCodeBlocks();

        for (var i = 0; i < codeBlocks.length; ++i) {
            operation(codeBlocks[i], "wrapped");
        }
    }

    function toggleLineWrapping() {
        var checkBox = getCheckBox();

        if (checkBox.checked) {
            forAllCodeBlocks(addClass);
        } else {
            forAllCodeBlocks(removeClass);
        }
    }

    function initControls() {
        if (findCodeBlocks().length > 0) {
            var checkBox = getCheckBox();
            var label = getLabelForCheckBox();

            checkBox.onclick = toggleLineWrapping;
            checkBox.checked = false;

            removeClass(label, "hidden");
         }
    }

    function switchTab() {
        var id = this.id.substr(1);

        for (var i = 0; i < tabs.tabs.length; i++) {
            if (tabs.tabs[i].id === id) {
                tabs.select(i);
                break;
            }
        }

        return false;
    }

    function select(i) {
        this.deselectAll();

        changeElementClass(this.tabs[i], "tab selected");
        changeElementClass(this.headers[i], "selected");

        while (this.headers[i].firstChild) {
            this.headers[i].removeChild(this.headers[i].firstChild);
        }

        var h2 = document.createElement("H2");

        h2.appendChild(document.createTextNode(this.titles[i]));
        this.headers[i].appendChild(h2);
    }

    function deselectAll() {
        for (var i = 0; i < this.tabs.length; i++) {
            changeElementClass(this.tabs[i], "tab deselected");
            changeElementClass(this.headers[i], "deselected");

            while (this.headers[i].firstChild) {
                this.headers[i].removeChild(this.headers[i].firstChild);
            }

            var a = document.createElement("A");

            a.setAttribute("id", "ltab" + i);
            a.setAttribute("href", "#tab" + i);
            a.onclick = switchTab;
            a.appendChild(document.createTextNode(this.titles[i]));

            this.headers[i].appendChild(a);
        }
    }

    function findTabs(container) {
        return findChildElements(container, "DIV", "tab");
    }

    function findHeaders(container) {
        var owner = findChildElements(container, "UL", "tabLinks");
        return findChildElements(owner[0], "LI", null);
    }

    function findTitles(tabs) {
        var titles = [];

        for (var i = 0; i < tabs.length; i++) {
            var tab = tabs[i];
            var header = findChildElements(tab, "H2", null)[0];

            header.parentNode.removeChild(header);

            if (header.innerText) {
                titles.push(header.innerText);
            } else {
                titles.push(header.textContent);
            }
        }

        return titles;
    }

    function findChildElements(container, name, targetClass) {
        var elements = [];
        var children = container.childNodes;

        for (var i = 0; i < children.length; i++) {
            var child = children.item(i);

            if (child.nodeType === 1 && child.nodeName === name) {
                if (targetClass && child.className.indexOf(targetClass) < 0) {
                    continue;
                }

                elements.push(child);
            }
        }

        return elements;
    }

    // Entry point.

    window.onload = function() {
        initTabs();
        initControls();
    };
} (window, window.document));
""")
        writer.flush()
        writer.close()
    }

    private static void writeCss(Path reportDir) {
        final indexFile = reportDir.resolve("style.css")
        final Charset charset = Charset.defaultCharset()
        Writer writer = Files.newBufferedWriter(indexFile, charset)

        writer.write("""
body {
    margin: 0;
    padding: 0;
    font-family: sans-serif;
    font-size: 12pt;
}

body, a, a:visited {
    color: #303030;
}

#content {
    padding-left: 50px;
    padding-right: 50px;
    padding-top: 30px;
    padding-bottom: 30px;
}

#content h1 {
    font-size: 160%;
    margin-bottom: 10px;
}

#footer {
    margin-top: 100px;
    font-size: 80%;
    white-space: nowrap;
}

#footer, #footer a {
    color: #a0a0a0;
}

#line-wrapping-toggle {
    vertical-align: middle;
}

#label-for-line-wrapping-toggle {
    vertical-align: middle;
}

ul {
    margin-left: 0;
}

h1, h2, h3 {
    white-space: nowrap;
}

h2 {
    font-size: 120%;
}

ul.tabLinks {
    padding-left: 0;
    padding-top: 10px;
    padding-bottom: 10px;
    overflow: auto;
    min-width: 800px;
    width: auto !important;
    width: 800px;
}

ul.tabLinks li {
    float: left;
    height: 100%;
    list-style: none;
    padding-left: 10px;
    padding-right: 10px;
    padding-top: 5px;
    padding-bottom: 5px;
    margin-bottom: 0;
    -moz-border-radius: 7px;
    border-radius: 7px;
    margin-right: 25px;
    border: solid 1px #d4d4d4;
    background-color: #f0f0f0;
}

ul.tabLinks li:hover {
    background-color: #fafafa;
}

ul.tabLinks li.selected {
    background-color: #c5f0f5;
    border-color: #c5f0f5;
}

ul.tabLinks a {
    font-size: 120%;
    display: block;
    outline: none;
    text-decoration: none;
    margin: 0;
    padding: 0;
}

ul.tabLinks li h2 {
    margin: 0;
    padding: 0;
}

div.tab {
}

div.selected {
    display: block;
}

div.deselected {
    display: none;
}

div.tab table {
    min-width: 350px;
    width: auto !important;
    width: 350px;
    border-collapse: collapse;
}

div.tab th, div.tab table {
    border-bottom: solid #d0d0d0 1px;
}

div.tab th {
    text-align: left;
    white-space: nowrap;
    padding-left: 6em;
}

div.tab th:first-child {
    padding-left: 0;
}

div.tab td {
    white-space: nowrap;
    padding-left: 6em;
    padding-top: 5px;
    padding-bottom: 5px;
}

div.tab td:first-child {
    padding-left: 0;
}

div.tab td.numeric, div.tab th.numeric {
    text-align: right;
}

span.code {
    display: inline-block;
    margin-top: 0em;
    margin-bottom: 1em;
}

span.code pre {
    font-size: 11pt;
    padding-top: 10px;
    padding-bottom: 10px;
    padding-left: 10px;
    padding-right: 10px;
    margin: 0;
    background-color: #f7f7f7;
    border: solid 1px #d0d0d0;
    min-width: 700px;
    width: auto !important;
    width: 700px;
}

span.wrapped pre {
    word-wrap: break-word;
    white-space: pre-wrap;
    word-break: break-all;
}

label.hidden {
    display: none;
}

#summary {
    margin-top: 30px;
    margin-bottom: 40px;
}

#summary table {
    border-collapse: collapse;
}

#summary td {
    vertical-align: top;
}

.breadcrumbs, .breadcrumbs a {
    color: #606060;
}

.infoBox {
    width: 110px;
    padding-top: 15px;
    padding-bottom: 15px;
    text-align: center;
}

.infoBox p {
    margin: 0;
}

.counter, .percent {
    font-size: 120%;
    font-weight: bold;
    margin-bottom: 8px;
}

#duration {
    width: 125px;
}

#successRate, .summaryGroup {
    border: solid 2px #d0d0d0;
    -moz-border-radius: 10px;
    border-radius: 10px;
}

#successRate {
    width: 140px;
    margin-left: 35px;
}

#successRate .percent {
    font-size: 180%;
}

.success, .success a {
    color: #008000;
}

div.success, #successRate.success {
    background-color: #bbd9bb;
    border-color: #008000;
}

.failures, .failures a {
    color: #b60808;
}

.skipped, .skipped a {
    color: #c09853;
}

div.failures, #successRate.failures {
    background-color: #ecdada;
    border-color: #b60808;
}

ul.linkList {
    padding-left: 0;
}

ul.linkList li {
    list-style: none;
    margin-bottom: 5px;
}

div.test h3 {
    display: inline;
}

span.code {
    display: block;
    margin-top: 10px;
    margin-bottom: 30px;
}

a.workdir {
    display: inline;
    margin-left: 5px;
}

""")

        writer.flush()
        writer.close()
    }


}
