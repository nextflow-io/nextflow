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


/**
 * Split a string into literal tokens i.e. sequence of chars that does not contain a blank
 * or wrapper by a quote or double quote
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
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