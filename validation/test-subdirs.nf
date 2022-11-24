workflow {
  foo | bar
}

process foo {
    output:
    path 'gsfolder/'
    path 'gsfolder2'
    path 'gsfolder3'
    path 'gsfolder4'
    path 'test5.txt'
    path 'test6.txt'
    path 'test7.txt'
    path 'gsfolder5/sub'

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

    mkdir -p gsfolder5/sub
    touch gsfolder5/sub/test8.txt
    """
}

process bar {
    input:
    path test_folder
    path test_folder2
    path 'test-folder3'
    path 'test-folder4/*'
    path test5
    path 'this-is-test-6.txt'
    path 'test7/foo/*'
    path 'test8/*'

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
