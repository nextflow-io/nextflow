
import groovyx.gpars.dataflow.DataflowVariable

final DataflowVariable variable = new DataflowVariable()
final DataflowVariable result = new DataflowVariable()

variable.then {it * 2} then {it + 1} then {result << it}
variable << 4

println result.val

