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

package nextflow.util

import groovy.transform.CompileStatic

/**
 * A service loader inspired to {@link ServiceLoader}
 * that allows to load only service names or classes without instantiating
 * the actual instances
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ServiceDiscover<S> {

    private static final String PREFIX = "META-INF/services/";

    private Class<S> service

    private ClassLoader classLoader

    static <S> List<Class<S>> load(Class<S> service) {
        new ServiceDiscover<S>(service).load()
    }

    ServiceDiscover(Class<S> service) {
        assert service
        this.service = service
        this.classLoader = Thread.currentThread().getContextClassLoader()
    }

    ServiceDiscover(Class<S> service, ClassLoader classLoader) {
        assert service
        assert classLoader
        this.service = service
        this.classLoader = classLoader
    }

    private int parseLine(Class service, URL u, BufferedReader r, int lc, List<String> names)
            throws IOException, ServiceConfigurationError
    {
        String ln = r.readLine();
        if (ln == null) {
            return -1
        }

        int ci = ln.indexOf('#')
        if (ci >= 0)
            ln = ln.substring(0, ci)

        ln = ln.trim()
        int n = ln.length()

        if (n != 0) {
            if ((ln.indexOf(' ') >= 0) || (ln.indexOf('\t') >= 0))
                fail(service, u, lc, "Illegal configuration-file syntax")

            int cp = ln.codePointAt(0)
            if (!Character.isJavaIdentifierStart(cp))
                fail(service, u, lc, "Illegal provider-class name: " + ln)

            for (int i = Character.charCount(cp); i < n; i += Character.charCount(cp)) {
                cp = ln.codePointAt(i)
                if (!Character.isJavaIdentifierPart(cp) && (cp != '.'))
                    fail(service, u, lc, "Illegal provider-class name: " + ln)
            }

            if ( !names.contains(ln) )
                names.add(ln)
        }

        return lc + 1
    }


    private List<String> parse(Class service, URL u) throws ServiceConfigurationError
    {
        InputStream stream = null;
        BufferedReader reader = null;
        ArrayList<String> names = new ArrayList<>();
        try {
            stream = u.openStream();
            reader = new BufferedReader(new InputStreamReader(stream, "utf-8"));
            int lc = 1;
            while ((lc = parseLine(service, u, reader, lc, names)) >= 0) { /*empty*/ }
        }
        catch (IOException x) {
            fail(service, "Error reading configuration file", x);
        }
        finally {
            try {
                if (reader != null) reader.close();
                if (stream != null) stream.close();
            }
            catch (IOException y) {
                fail(service, "Error closing configuration file", y);
            }
        }
        return names
    }

    private static void fail(Class service, String msg, Throwable cause)
            throws ServiceConfigurationError
    {
        throw new ServiceConfigurationError(service.getName() + ": " + msg,
                cause);
    }

    private static void fail(Class service, String msg)
            throws ServiceConfigurationError
    {
        throw new ServiceConfigurationError(service.getName() + ": " + msg);
    }

    private static void fail(Class service, URL u, int line, String msg)
            throws ServiceConfigurationError
    {
        fail(service, u.toString() + ":" + line + ": " + msg);
    }

    private Class<S> classForName(String clazzName) {

        Class<?> result = null
        try {
            result = Class.forName(clazzName, false, classLoader);
        }
        catch (ClassNotFoundException x) {
            fail(service, "Provider $clazzName not found");
        }
        if (!service.isAssignableFrom(result)) {
            fail(service, "Provider $clazzName not a subtype");
        }
        try {
            return (Class<S>)result;
        }
        catch (Throwable x) {
            fail(service, "Provider $clazzName could not be instantiated: " + x, x);
        }
        throw new Error()
    }

    List<String> names() {

        Enumeration<URL> configs = null

        try {
            configs = classLoader.getResources(PREFIX + service.getName())

        } catch (IOException x) {
            fail(service, "Error locating configuration files", x);
        }

        List<String> result = []
        for( URL url : configs ) {
            for( String name : parse(service, url) ) {
                result << name
            }
        }
        return result
    }

    List<Class<S>> load() {
        names().collect{ classForName(it) }
    }


}
