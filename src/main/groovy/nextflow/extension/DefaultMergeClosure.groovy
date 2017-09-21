package nextflow.extension
import groovyx.gpars.dataflow.operator.DataflowProcessor
/**
 * Implement the default `merge` operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultMergeClosure extends Closure {

    private int numOfParams

    DefaultMergeClosure(int n) {
        super(null, null);
        numOfParams = n
    }

    @Override
    public int getMaximumNumberOfParameters() {
         numOfParams
    }

    @Override
    public Class[] getParameterTypes() {
        Collections.nCopies(numOfParams, Object)
    }

    @Override
    public void setDelegate(final Object delegate) {
        super.setDelegate(delegate);
    }

    @Override
    public void setResolveStrategy(final int resolveStrategy) {
        super.setResolveStrategy(resolveStrategy);
    }

    @Override
    public Object call(final Object arguments) {
        throw new UnsupportedOperationException()
    }

    @Override
    public Object call(final Object... args) {
        final result = []
        for( int i=0; i<args.size(); i++ )
            DataflowHelper.addToList(result, args[i])
        ((DataflowProcessor) getDelegate()).bindAllOutputsAtomically(result);
        return result;
    }

    @Override
    public Object call() {
        throw new UnsupportedOperationException()
    }

}
