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

package nextflow.util
import java.nio.file.Files
import java.nio.file.Path
import java.util.regex.Matcher
import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Read a INI file
 *
 * See http://stackoverflow.com/a/15638381/395921
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class IniFile {

    private Pattern fSection = Pattern.compile( "\\s*\\[([^]]*)\\]\\s*" );
    private Pattern fKeyValue = Pattern.compile( "\\s*([^=]*)=(.*)" );
    private Map<String, Map<String, String>> fEntries = new HashMap<>();

    private Path fFile;

    IniFile() {}

    IniFile( Path path ) throws IOException {
        load( path );
    }

    IniFile(String path)  {
        load(path)
    }

    IniFile load( String path ) {
        assert path
        load(path as Path)
    }

    IniFile load( File file ) {
        assert file
        load(file.toPath())
    }

    IniFile load( Path path, cs = null ) {
        assert path

        this.fFile = path
        if( Files.exists(path) ) {
            final charset = CharsetHelper.getCharset(cs)
            final reader = Files.newBufferedReader(path, charset)
            try {
                load(reader)
            }
            finally{
                reader.close()
            }
        }

        return this
    }

    IniFile load( Reader br ) {
        assert br

        String line;
        String section = null;
        while(( line = br.readLine()) != null ) {
            Matcher m = fSection.matcher( line );
            if( m.matches()) {
                section = m.group( 1 ).trim();
            }
            else if( section != null ) {
                m = fKeyValue.matcher( line );
                if( m.matches()) {
                    String key   = m.group( 1 ).trim();
                    String value = m.group( 2 ).trim();
                    Map< String, String > kv = fEntries.get( section );
                    if( kv == null ) {
                        fEntries.put( section, kv = new HashMap<>());
                    }
                    kv.put( key, value );
                }
            }
        }

        return this
    }

    String getString( String section, String key, String defValue = null ) {
        Map< String, String > kv = fEntries.get( section );
        if( kv == null ) {
            return defValue;
        }
        return kv.get(key) ?: defValue
    }

    int getInt( String section, String key, int defValue = 0) {
        Map< String, String > kv = fEntries.get( section );
        if( kv == null ) {
            return defValue;
        }

        kv.containsKey(key) ? Integer.parseInt( kv.get( key )) : defValue;
    }

    float getFloat( String section, String key, float defValue = 0 ) {
        Map< String, String > kv = fEntries.get( section );
        if( kv == null ) {
            return defValue;
        }

        kv.containsKey(key) ?  Float.parseFloat( kv.get( key )) : defValue
    }

    double getDouble( String section, String key, double defValue = 0 ) {
        Map< String, String > kv = fEntries.get( section );
        if( kv == null ) {
            return defValue;
        }

        kv.containsKey(key) ? Double.parseDouble( kv.get( key )) : defValue
    }

    boolean getBool( String section, String key, boolean defValue = false ) {
        Map< String, String > kv = fEntries.get( section );
        if( kv == null ) {
            return defValue;
        }

        kv.containsKey(key) ? Boolean.parseBoolean( kv.get( key )) : defValue
    }

    Map<String,String> section(String section) {
        def result = fEntries.get(section)
        return result != null ? Collections.unmodifiableMap(result) : Collections.<String,String>emptyMap()
    }

    def propertyMissing(String name) {
        if( fEntries.containsKey(name) )
            return section(name)

        throw new MissingPropertyException(name,IniFile)
    }

    def getFile() {
        return fFile
    }
}