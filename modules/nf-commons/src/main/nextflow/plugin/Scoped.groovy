/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.plugin

import java.lang.annotation.ElementType
import java.lang.annotation.Retention
import java.lang.annotation.RetentionPolicy
import java.lang.annotation.Target

/**
 * Allow the definition of plugin priority order, smaller value is an higher priority
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Deprecated
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
@interface Scoped {
    String value()
    int priority() default 0
}
