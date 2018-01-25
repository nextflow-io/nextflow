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

package nextflow.trace

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * Model a process summary data used to render box-plots in the execution HTML report
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class ReportSummary {

    private static Map<String,Closure<Double>> mappers = [:]

    static {
        mappers.cpu = { TraceRecord record -> record.get('%cpu') as Double }
        mappers.mem = { TraceRecord record -> record.get('vmem') as Double }
        mappers.time = { TraceRecord record -> record.get('realtime') as Double }
        mappers.reads = { TraceRecord record -> record.get('read_bytes') as Double }
        mappers.writes = { TraceRecord record -> record.get('write_bytes') as Double}
    }

    /**
     * Hold the summary for each series ie. cpu, memory, time, disk reads, disk writes
     */
    private Map<String, Summary> series

    /**
     * List of series names
     */
    private List<String> names = Collections.unmodifiableList(mappers.keySet() as List)

    ReportSummary() {
        this.series = new LinkedHashMap<>()
        for( String name : names ) {
            series[name] = new Summary(mappers[name])
        }
    }

    /**
     * @return The list of series names eg `cpu`, `time`, etc
     */
    List<String> getNames() { names }

    /**
     * Add a {@link TraceRecord} for which compute the final summary
     *
     * @param record A {@link TraceRecord} representing a computed task
     */
    void add( TraceRecord record ) {
        for( int i=0; i< names.size(); i++ ) {
            final name = names[i]
            series[name].add(record)
        }
    }

    /**
     * Compute the summary stats for the collected tasks
     *
     * @param
     *      name A series name eg {@code cpu}, {@code time}, etc
     * @return
     *      A {@link Map} holding the summary containing the following stats:
     *      - min: minimal value
     *      - q1: first quartile
     *      - q2: second quartile ie. median
     *      - q3: third quartile
     *      - max: maximum value
     *      - mean: average value
     *      - minLabel: label fot the task reporting the min value
     *      - q1Label: label fot the task reporting the q1 value
     *      - q2Label: label fot the task reporting the q2 value
     *      - q3Label: label fot the task reporting the q3 value
     *      - maxLabel: label fot the task reporting the max value
     */
    Map<String,Double> compute(String name) {
        if( !names.contains(name) )
            throw new IllegalArgumentException("Invalid status status field name: $name -- it must be one of the following: ${names.join(',')}")

        return series[name].compute()
    }

    /**
     * Model and compute a series summary
     */
    @Slf4j
    @CompileStatic
    static class Summary {

        private List<TraceRecord> tasks

        private Closure<Double> metric

        private BigDecimal total = 0

        private int count

        private Integer digits = 2

        @PackageScope Double min

        @PackageScope Double max

        @PackageScope String minLabel

        @PackageScope String maxLabel

        @PackageScope String q1Label

        @PackageScope String q2Label

        @PackageScope String q3Label

        Summary( Closure<Double> metric )  {
            this.metric = metric
            this.tasks = []
            this.min
        }

        private void label( def r, int q) {
            if( r instanceof TraceRecord ) {
                switch(q) {
                    case 0: minLabel = r.get('name'); break
                    case 25: q1Label = r.get('name'); break
                    case 50: q2Label = r.get('name'); break
                    case 75: q3Label = r.get('name'); break
                    case 100: maxLabel = r.get('name'); break
                    default:
                        log.debug "Invalid summary stats quantile: $q"
                }
            }
        }

        private double round( double value ) {
            if( digits == null )
                return value
            final x = Math.pow(10,digits)
            Math.round( value * x ) / x
        }

        void add( TraceRecord record ) {
            final value = metric(record)
            if( value == null )
                return
            count++
            tasks.add(record)
            total += value
            if( min==null || value<min ) {
                min = value
                minLabel = record.get('name')
            }
            if( max==null || value>max ) {
                max = value
                maxLabel = record.get('name')
            }
        }

        /**
         * Compute the stats for the collected tasks
         *
         * @return
         *      A {@link Map} holding the summary containing the following stats:
         *      - min: minimal value
         *      - q1: first quartile
         *      - q2: second quartile ie. median
         *      - q3: third quartile
         *      - max: maximum value
         *      - mean: average value
         *      - minLabel: label fot the task reporting the min value
         *      - q1Label: label fot the task reporting the q1 value
         *      - q2Label: label fot the task reporting the q2 value
         *      - q3Label: label fot the task reporting the q3 value
         *      - maxLabel: label fot the task reporting the max value
         */
        Map<String,Double> compute() {
            if( count==0 )
                return null

            final result = [:]
            final sorted = tasks.sort( false, { TraceRecord record -> metric(record) } )

            result.mean = round(total / count as double)
            result.min = round(quantile(sorted, 0))
            result.q1 = round(quantile(sorted, 25))
            result.q2 = round(quantile(sorted, 50))
            result.q3 = round(quantile(sorted, 75))
            result.max = round(quantile(sorted, 100))
            result.minLabel = minLabel
            result.maxLabel = maxLabel
            result.q1Label = q1Label
            result.q2Label = q2Label
            result.q3Label = q3Label
            return result
        }

        /**
         * Calculate the q-th quantile in the same way as the R `summary` function does
         * See https://www.r-bloggers.com/exploratory-data-analysis-the-5-number-summary-two-different-methods-in-r/
         *
         * @param items A list of numbers
         * @param q The q-th quantile. It must be a number between 0 and 100
         * @return The q-th quantile value
         */
        protected double quantile(List items, int q) {
            assert items, 'Argument items cannot be empty'
            assert q>=0 && q<=100
            if( q==0 ) {
                final X = items[0]
                label(X,q)
                return metric(X)
            }
            if( q==100 ) {
                final X = items[items.size()-1]
                label(X,q)
                return metric(X)
            }

            int n = items.size()
            double j = (n-1) * q / 100
            if( j == Math.floor(j) ) {
                final item = items[(int)j]
                label(item,q)
                return metric(item)
            }
            else {
                int i = (int)Math.floor(j); int k = (int)Math.ceil(j)
                label(items[i],q)
                final Xi = metric(items[i])
                final Xk = metric(items[k])
                return Xi + (j-i) * (Xk-Xi)
            }
        }

        protected double[] quantiles(List items) {
            def result = new double[5]
            result[0] = quantile(items,0)
            result[1] = quantile(items,25)
            result[2] = quantile(items,50)
            result[3] = quantile(items,75)
            result[4] = quantile(items,100)
            return result
        }
    }

}
