process foo {
    output:
    path 'gsfolder/' into ch1
    path 'gsfolder2' into ch2
    path 'gsfolder3' into ch3
    path 'gsfolder4' into ch4
    path 'test5.txt' into ch5
    path 'test6.txt' into ch6
    path 'test7.txt' into ch7

    """
    mkdir -p gsfolder/sub
    touch gsfolder/test1.txt
    touch gsfolder/sub/test1.txt

    mkdir gsfolder2
    touch gsfolder2/test2.txt

    mkdir gsfolder3
    touch gsfolder3/test3.txt

    mkdir gsfolder4
    touch gsfolder4/test4.txt

    echo 'Hello 5'>>test5.txt
    echo 'Hello 6'>>test6.txt
    echo 'Hello 7'>>test7.txt

    """
}

process bar {
    input:
    path test_folder from ch1
    path test_folder2 from ch2
    path 'test-folder3' from ch3
    path 'test-folder4/*' from ch4
    path test5 from ch5
    path 'this-is-test-6.txt' from ch6
    path 'test7/foo/*' from ch7

    """
    set -x
    [[ -f $test_folder/test1.txt ]] || false
    [[ -f $test_folder/sub/test1.txt ]] || false
    [[ -f $test_folder2/test2.txt ]] || false
    [[ -f test-folder3/test3.txt ]] || false
    [[ -f test-folder4/gsfolder4/test4.txt ]] || false
    [[ \$(cat $test5) = 'Hello 5' ]] || false
    [[ \$(cat this-is-test-6.txt) = 'Hello 6' ]] || false
    [[ \$(cat test7/foo/test7.txt) = 'Hello 7' ]] || false
    """
}
