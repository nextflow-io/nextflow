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

package nextflow.ui

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
import groovy.transform.EqualsAndHashCode

@EqualsAndHashCode(includes='value')
class TextLabel {

    enum Align { RIGHT, LEFT }

    int width

    Align align = Align.LEFT

    boolean active

    def value

    def List<LabelDecorator> decorators = []

    def Integer max

    /**
     * Create a label with the provided value
     */
    TextLabel( def value ) {
        this.value = value
        this.align = ( value instanceof Number || value?.toString()?.isNumber() ) ? Align.RIGHT : Align.LEFT
    }

    def TextLabel width( int num ) { this.width = num; this }

    def TextLabel left() { this.align = Align.LEFT; this }

    def TextLabel right() { this.align = Align.RIGHT; this }

    def TextLabel number() { this.align = Align.RIGHT; this }

    def TextLabel max( int value ) { this.max = value; this }

    /**
     * Switch OFF the decorators rendering
     */
    def TextLabel switchOff() { active = false; this }

    /**
     * Turn ON the decorators rendering
     */
    def TextLabel switchOn() {
        if ( !decorators ) {
            decorators << AnsiStyle.style().negative()
        }

        active = true;
        return this
    }


    def TextLabel leftShift( LabelDecorator deco ) {
        this.decorators.add(deco)
        return this
    }

    def TextLabel add( LabelDecorator deco ) {
        this.decorators.add(deco)
        return this
    }

    /**
     * @return Renders the string applying the provided decorators
     */
    def String toString() {

        String result = value != null ? value.toString() : '-'

        if( width && align == Align.LEFT ) {
            result = result.padRight(width)
        }
        else if ( width && align == Align.RIGHT ) {
            result = result.padLeft(width)
        }
        else {
            result
        }

        // apply the max rule
        if ( max != null ) {
            result = applyMax(result)
        }

        // if not active return as it is
        if( !active )  {
            return result
        }

        decorators.each {
            result = it.apply(this,result)
        }

        result
    }

    String applyMax( String str ) {
        assert str

        if( str.length()<=max) {
            return str
        }

        def cut
        if( align == Align.LEFT ) {
            cut = str.substring(max)
            str = str.substring(0,max)
        }
        else {
            int p = str.size()-max
            cut = str.substring( 0, p )
            str = str.substring( p, str.size() )
        }

        if( str.size()>3 && !cut.isAllWhitespace() ) {
            if ( align == Align.LEFT ) {
                str = str[0..-3] + '..'
            }
            else if( align == Align.RIGHT  ) {
                str = '..' + str[2..-1]
            }
        }


        return str

    }

    /**
     * Create a new {@code TextLabel} object with the specified value
     */
    static TextLabel of( def value ) {
        new TextLabel( value )
    }

}