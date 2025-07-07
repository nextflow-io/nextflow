/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.types.shim;

import nextflow.script.dsl.Description;

@Description("""
    A string is an immutable sequence of characters.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#string)
""")
@ShimType(java.lang.String.class)
public interface String {

    @Description("""
        Returns `true` if the string ends with the given suffix.
    """)
    boolean endsWith(String suffix);

    @Description("""
        Execute the string as a command. Returns a [Process](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Process.html) which provides the exit status and standard input/output/error of the executed command.
    """)
    Process execute();

    @Description("""
        Returns the index within the string of the first occurrence of the given substring. Returns -1 if the string does not contain the substring.
    """)
    int indexOf(String str);

    @Description("""
        Returns the index within the string of the first occurrence of the given substring, starting the search at the given index. Returns -1 if the string does not contain the substring.
    """)
    int indexOf(String str, int fromIndex);

    @Description("""
        Returns `true` if the string is empty or contains only whitespace characters.
    """)
    boolean isBlank();

    @Description("""
        Returns `true` if the string is empty (i.e. `length()` is 0).
    """)
    boolean isEmpty();

    @Description("""
        Returns `true` if the string can be parsed as a floating-point number.
    """)
    boolean isFloat();

    @Description("""
        Returns `true` if the string can be parsed as an integer.
    """)
    boolean isInteger();

    @Description("""
        Returns the index within the string of the last occurrence of the given substring. Returns -1 if the string does not contain the substring.
    """)
    int lastIndexOf(String str);

    @Description("""
        Returns the index within the string of the last occurrence of the given substring, searching backwards starting at the given index. Returns -1 if the string does not contain the substring.
    """)
    int lastIndexOf(String str, int fromIndex);

    @Description("""
        Returns the length of the string.
    """)
    int length();

    @Description("""
        Returns the MD5 checksum of the string.
    """)
    String md5();

    @Description("""
        Returns a new string in which each occurrence of the target string is replaced with the given replacement string.
    """)
    String replace(String target, String replacement);

    @Description("""
        Returns a new string in which each occurrence of the given regular expression is replaced with the given replacement string.
    """)
    String replaceAll(String regex, String replacement);

    @Description("""
        Returns a new string in which the first occurrence of the given regular expression is replaced with the given replacement string.
    """)
    String replaceFirst(String regex, String replacement);

    @Description("""
        Returns the SHA-256 checksum of the string.
    """)
    String sha256();

    @Description("""
        Returns `true` if the string ends with the given prefix.
    """)
    boolean startsWith(String prefix);

    @Description("""
        Returns a copy of the string with all leading and trailing whitespace removed.
    """)
    String strip();

    @Description("""
        Returns a copy of the string with leading spaces on each line removed. The number of spaces to remove is determined by the line with the least number of leading spaces, excluding lines with only whitespace.
    """)
    String stripIndent();

    @Description("""
        Returns a copy of the string with all leading whitespace removed.
    """)
    String stripLeading();

    @Description("""
        Returns a copy of the string with all trailing whitespace removed.
    """)
    String stripTrailing();

    @Description("""
        Returns a substring of this string.
    """)
    String substring(int beginIndex);

    @Description("""
        Returns a substring of this string.
    """)
    String substring(int beginIndex, int endIndex);

    @Description("""
        Returns `true` if the trimmed string is "true", "y", or "1" (ignoring case).
    """)
    Boolean toBoolean();

    @Description("""
        Parses the string into a floating-point number.
    """)
    Float toFloat();

    @Description("""
        Parses the string into an integer.
    """)
    Integer toInteger();

    @Description("""
        Returns a copy of this string with all characters converted to lower case.
    """)
    String toLowerCase();

    @Description("""
        Returns a copy of this string with all characters converted to upper case.
    """)
    String toUpperCase();

    @Description("""
        Splits the string into a list of substrings using the given delimiters. Each character in the delimiter string is treated as a separate delimiter.
    """)
    List<String> tokenize(String delimiters);

}
