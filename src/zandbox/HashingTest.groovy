import embed.com.google.common.hash.HashCodes
import nextflow.util.CacheHelper

/*
 * Copyright (c) 2012, the authors.
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

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

def hash = CacheHelper.hasher(0).hash()
println hash.toString()
println hash.asBytes()
println HashCodes.fromBytes(hash.asBytes()).toString()


assert CacheHelper.hasher(['abc',123]).hash() == CacheHelper.hasher(['abc',123]).hash()
assert CacheHelper.hasher(['abc',123]).hash() != CacheHelper.hasher([123,'abc']).hash()