process absoluteExternalPath {

    output:
    path "absoluteExternalPath.txt" into absoluteExternalPathOut

    """
    ln -s $baseDir/data/p1.fa absoluteExternalPath.txt
    """
}
absoluteExternalPathOut.view{ assert it.text == ("$baseDir/data/p1.fa" as Path).text; "absoluteExternalPath TRUE" }

process relativeExternalPath {

    output:
    path "relativeExternalPath.txt" into relativeExternalPathOut

    """
    link=`realpath -s --relative-to="\$PWD" "$baseDir/data/p1.fa"`
    echo current: \$PWD
    echo rellink: \$link
    ln -s \$link relativeExternalPath.txt
    """
}
relativeExternalPathOut.view{ assert it.text == ("$baseDir/data/p1.fa" as Path).text; "relativeExternalPath TRUE" }

process absoluteInternalPath {

    output:
    path "absoluteInternalPath.txt" into absoluteInternalPathOut

    """
    printf "absoluteInternalPathText" > absoluteInternalPathFile.txt
    ln -s \$PWD/absoluteInternalPathFile.txt absoluteInternalPath.txt
    """
}
absoluteInternalPathOut.view{ assert it.text == "absoluteInternalPathText"; "absoluteInternalPath TRUE" }

process absoluteInternalPathInDir {

    output:
    path "absoluteInternalPathInDir.txt" into absoluteInternalPathOutInDir

    """
    mkdir a
    printf "absoluteInternalPathInDirText" > a/absoluteInternalPathFileInDir.txt
    ln -s \$PWD/a/absoluteInternalPathFileInDir.txt absoluteInternalPathInDir.txt
    """
}
absoluteInternalPathOutInDir.view{ assert it.text == "absoluteInternalPathInDirText"; "absoluteInternalPathInDir TRUE" }

process relativeInternalPath {

    output:
    path "relativeInternalPath.txt" into relativeInternalPathOut

    """
    printf "relativeInternalPathText" > relativeInternalPathFile.txt
    ln -s relativeInternalPathFile.txt relativeInternalPath.txt
    """
}
relativeInternalPathOut.view{ assert it.text == "relativeInternalPathText"; "relativeInternalPath TRUE" }

process relativeInternalPathInDir {

    output:
    path "relativeInternalPathInDir.txt" into relativeInternalPathOutInDir

    """
    mkdir a
    printf "relativeInternalPathInDirText" > a/relativeInternalPathFileInDir.txt
    ln -s a/relativeInternalPathFileInDir.txt relativeInternalPathInDir.txt
    """
}
relativeInternalPathOutInDir.view{ assert it.text == "relativeInternalPathInDirText"; "relativeInternalPathInDir TRUE" }

process relativeInternalPathBothInDir {

    output:
    path "b/relativeInternalPathBothInDir.txt" into relativeInternalPathOutBothInDir

    """
    mkdir a b
    printf "relativeInternalPathBothInDirText" > a/relativeInternalPathFileBothInDir.txt
    ln -s ../a/relativeInternalPathFileBothInDir.txt b/relativeInternalPathBothInDir.txt
    """
}
relativeInternalPathOutBothInDir.view{ assert it.text == "relativeInternalPathBothInDirText"; "relativeInternalPathBothInDir TRUE" }

process relativeInternalPathBothInDirButMoved {

    output:
    path "a/"
    path "b/relativeInternalPathBothInDirButMoved.txt" into relativeInternalPathOutBothInDirButMoved

    """
    mkdir a b
    printf "relativeInternalPathFileBothInDirButMovedText" > a/relativeInternalPathFileBothInDirButMoved.txt
    ln -s ../a/relativeInternalPathFileBothInDirButMoved.txt b/relativeInternalPathBothInDirButMoved.txt
    """
}
relativeInternalPathOutBothInDirButMoved.view{ assert it.text == "relativeInternalPathFileBothInDirButMovedText"; "relativeInternalPathFileBothInDirButMoved TRUE" }

process linkInFolder {

    output:
    path "a" into linkInFolderOut

    """
    mkdir a
    mkdir a/b
    printf "linkInFolderText" > linkInFolderFile.txt
    ln -s ../../linkInFolderFile.txt a/b/linkInFolder.txt
    """
}
linkInFolderOut.view { assert it.resolve("b/linkInFolder.txt").text == "linkInFolderText"; "linkInFolder TRUE" }