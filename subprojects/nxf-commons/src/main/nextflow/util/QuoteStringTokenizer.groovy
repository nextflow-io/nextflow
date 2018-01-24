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


/**
 * Split a string into literal tokens i.e. sequence of chars that does not contain a blank
 * or wrapper by a quote or double quote
 *
 *
 * @author Paolo Di Tommaso
 *
 */
public class QuoteStringTokenizer implements Iterator<String>, Iterable<String> {

    private List<Character> chars = Arrays.<Character>asList(' ' as Character);

    private List<String> tokens = new ArrayList<String>();

    private Iterator<String> itr;

    public QuoteStringTokenizer( String value ) {
        this(value, ' ' as Character);
    }

    public QuoteStringTokenizer( String value, char... separators ) {

        if( separators == null || separators.length==0 ) {
            chars = new ArrayList<Character>(3);
            chars .add(' ' as Character);
        }
        else {
            chars = new ArrayList<Character>(3);
            for( char ch : separators ) { chars.add(ch); }
        }

        chars.add("\"" as char);
        chars.add("'" as char);

        /* start parsing */
        parseNext(value != null ? value.trim() : "");

    }

    void parseNext( String value )  {
        for( int i=0; i<value.length(); i++  ) {

            char ch=value.charAt(i);
            if( chars.contains(ch) ) {
                if( i>0 ) {
                    tokens.add(value.substring(0,i));
                }

                if( ch=='"' as char || ch=='\'' as char ) {
                    parseQuote(value.substring(i+1), ch);
                }
                else if( chars.contains(ch) ) {
                    parseNext(value.substring(i+1));
                }
                break;
            }
            // and of the string
            else if( i+1 == value.length() ) {
                tokens.add(value);
            }
        }
    }

    void parseQuote( String value, char delim ) {

        int p = value.indexOf(delim as String);
        if( p == -1 ) {
            tokens.add(value);
            return;
        }

        tokens.add( value.substring(0,p) );
        parseNext(value.substring(p+1));
    }

    public boolean hasNext() {
        return itr().hasNext();
    }

    public String next() {
        return itr().next();
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove not supported");
    }

    public String toString() {
        return tokens.toString();
    }

    public Iterator<String> iterator() {
        return itr();
    }

    /*
      * lazy iterator creator
      */
    private Iterator<String> itr() {
        if( itr == null ) {
            itr = tokens.iterator();
        }
        return itr;
    }


}