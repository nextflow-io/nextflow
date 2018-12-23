/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
 * Used to render text based table for UI purpose. Example:
 *
 * <pre>
 *     def table = new TableBuilder().head('col1').head('col2').head('col3', 10)
 *     table << x << y << z << table.closeRow()
 *     table << p << q << w << table.closeRow()
 *
 *     :
 *     println table.toString()
 *     </pre>
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TableBuilder {

    /**
     * The list defining the table header
     */
    List<TextLabel> headers = []

    List<TextLabel.Align> colsAlign = []

    /**
     * All the rows
     */
    List<List<TextLabel>> allRows = []

    List<Integer> dim = []

    List<Integer> maxColsWidth = []

    String cellSeparator = ' '

    String rowSeparator = '\n'

    List<TextLabel> currentRow = []

    /**
     * Defines a single column header, for example:
     * <pre>
     *     def table = new TableBuilder().head('col1').head('col2').head('col3', 10)
     *     table << x << y << z
     *     </pre>
     *
     * @param name The string value to be used as column header
     * @param maxWidth The column max width
     * @return The table itself, to enable method chaining
     */
    TableBuilder head( String name, int maxWidth = 0 ) {
        def label = new TextLabel(name)

        headers << label
        maxColsWidth << maxWidth
        colsAlign << null

        def widths = headers .collect { it?.toString()?.size() }
        trackWidths(widths)
        return this
    }

    TableBuilder head( String name, TextLabel.Align align ) {
        assert align
        head(name)
        colsAlign[-1] = align

        return this
    }

    /**
     * Use the specific string array as the table header definition
     * @param cols
     * @return
     */
    TableBuilder setHeaders( String... cols ) {
        assert cols != null
        setHeaders( cols.collect { new TextLabel(it) } as TextLabel[] )
    }

    TableBuilder setHeaders( TextLabel... cols ) {
        assert cols != null
        // copy the header
        this.headers = new ArrayList<>(cols as List<TextLabel>)
        // keep track of the columns width
        def widths = cols .collect { it?.toString()?.size() }
        trackWidths(widths)
        //return the object itself
        return this
    }


    TableBuilder setMaxColsWidth( int...colsWidth ) {
        assert colsWidth != null

        maxColsWidth = new ArrayList<>(colsWidth as List<Integer>)

        return this
    }

    TableBuilder append( Object... values ) {
        append( values as List )
    }

    /**
     * Append the specified list of values as the next row in the table, the columns width
     * are adapted accordingly
     *
     * @param values
     * @return
     */
    TableBuilder append( List values ) {
        assert values != null

        def row = new ArrayList<TextLabel>(values.size())
        def len = new ArrayList<Integer>(values.size())
        values.each{  it ->
            row << ( it instanceof TextLabel ? it : new TextLabel(it) )
            len << it?.toString()?.size()
        }

        trackWidths(len)
        allRows << row

        this
    }

    /**
     * Defines the left-shift operator useful to build the table using the following syntax
     * <pre>
     *      def table = new TableBuilder()
     *      table << col1 << col2 << col3 << table.closeRow()
     *      :
     *     </pre>
     *
     *
     * @param value The value to be added in the table at the current row
     * @return The table instance itself
     */
    TableBuilder leftShift( def value ) {
        if( value == this ) {
            return this
        }

        if( value instanceof TextLabel ) {
            currentRow << value
        }
        else {
            currentRow << new TextLabel(value)
        }

        return this
    }

    /**
     * Close a row in the
     * @return
     */
    TableBuilder closeRow() {
        append(currentRow)
        currentRow = []
        return this
    }

    protected void trackWidths( List<Integer> newDim ) {
        def size = Math.min(dim.size(), newDim.size())

        for( int i=0; i<size; i++ ) {
            if ( dim[i] < newDim[i])  {
                dim[i] = newDim[i]
            }
        }

        if( newDim.size() > size ) {
            for( int i=size; i<newDim.size(); i++ ) {
                dim << newDim[i]
            }
        }
    }

    /**
     * @return Render the final table and return the string
     */
    String toString() {

        // check if there's some rows no closed
        if( currentRow ) {
            closeRow()
        }


        StringBuilder result = new StringBuilder()

        def count=0

        /*
         * render the header
         */
        if( headers ) {
            count++
            headers.eachWithIndex { TextLabel cell, int index ->
                renderCell( result, cell, index )
            }
        }


        /*
         * render the table
         */
        allRows.each{ List<TextLabel> row ->
            // render the 'rowSeparator' (only after the first row
            if( count++ && rowSeparator!=null )  { result.append(rowSeparator) }
            // render the row
            row.eachWithIndex { TextLabel cell, int index ->
                renderCell( result, cell, index )
            }
        }

        result.toString()
    }

    /**
     * Render a cell in the table
     *
     * @param result The {@code StringBuilder} collecting the result table text
     * @param cell The cell to the rendered
     * @param index The current index in the row of the cell to be rendered
     */
    private void renderCell( StringBuilder result, TextLabel cell, int index ) {

        // the 'cellSeparator' only after the first col
        if( index && cellSeparator != null ) result.append(cellSeparator)

        // set the max col width
        if( maxColsWidth && index<maxColsWidth.size() && maxColsWidth[index] ) {
            cell.max( maxColsWidth[index] )
        }

        // set the max
        cell.width( dim[index] )

        if( colsAlign[index] ) {
            cell.setAlign( colsAlign[index] )
        }

        // render the cell
        result.append( cell.toString() )

    }

}