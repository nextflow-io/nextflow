process absoluteExternalPath {

    output:
    path "absoluteExternalPath.txt"

    """
    ln -s $baseDir/data/p1.fa absoluteExternalPath.txt
    """
}

process relativeExternalPath {

    output:
    path "relativeExternalPath.txt"

    """
    link=`realpath -s --relative-to="\$PWD" "$baseDir/data/p1.fa"`
    echo current: \$PWD
    echo rellink: \$link
    ln -s \$link relativeExternalPath.txt
    """
}

process absoluteInternalPath {

    output:
    path "absoluteInternalPath.txt"

    """
    printf "absoluteInternalPathText" > absoluteInternalPathFile.txt
    ln -s \$PWD/absoluteInternalPathFile.txt absoluteInternalPath.txt
    """
}

process absoluteInternalPathInDir {

    output:
    path "absoluteInternalPathInDir.txt"

    """
    mkdir a
    printf "absoluteInternalPathInDirText" > a/absoluteInternalPathFileInDir.txt
    ln -s \$PWD/a/absoluteInternalPathFileInDir.txt absoluteInternalPathInDir.txt
    """
}

process relativeInternalPath {

    output:
    path "relativeInternalPath.txt"

    """
    printf "relativeInternalPathText" > relativeInternalPathFile.txt
    ln -s relativeInternalPathFile.txt relativeInternalPath.txt
    """
}

process relativeInternalPathInDir {

    output:
    path "relativeInternalPathInDir.txt"

    """
    mkdir a
    printf "relativeInternalPathInDirText" > a/relativeInternalPathFileInDir.txt
    ln -s a/relativeInternalPathFileInDir.txt relativeInternalPathInDir.txt
    """
}

process relativeInternalPathBothInDir {

    output:
    path "b/relativeInternalPathBothInDir.txt"

    """
    mkdir a b
    printf "relativeInternalPathBothInDirText" > a/relativeInternalPathFileBothInDir.txt
    ln -s ../a/relativeInternalPathFileBothInDir.txt b/relativeInternalPathBothInDir.txt
    """
}

process relativeInternalPathBothInDirButMoved {

    output:
    path "a/"
    path "b/relativeInternalPathBothInDirButMoved.txt"

    """
    mkdir a b
    printf "relativeInternalPathFileBothInDirButMovedText" > a/relativeInternalPathFileBothInDirButMoved.txt
    ln -s ../a/relativeInternalPathFileBothInDirButMoved.txt b/relativeInternalPathBothInDirButMoved.txt
    """
}

process linkInFolder {

    output:
    path "a"

    """
    mkdir a
    mkdir a/b
    printf "linkInFolderText" > linkInFolderFile.txt
    ln -s ../../linkInFolderFile.txt a/b/linkInFolder.txt
    """
}

workflow {
    absoluteExternalPath()
        .view{ assert it.text == ("$baseDir/data/p1.fa" as Path).text; "absoluteExternalPath TRUE" }
    relativeExternalPath()
        .view{ assert it.text == ("$baseDir/data/p1.fa" as Path).text; "relativeExternalPath TRUE" }
    absoluteInternalPath()
        .view{ assert it.text == "absoluteInternalPathText"; "absoluteInternalPath TRUE" }
    absoluteInternalPathInDir()
        .view{ assert it.text == "absoluteInternalPathInDirText"; "absoluteInternalPathInDir TRUE" }
    relativeInternalPath()
        .view{ assert it.text == "relativeInternalPathText"; "relativeInternalPath TRUE" }
    relativeInternalPathInDir()
        .view{ assert it.text == "relativeInternalPathInDirText"; "relativeInternalPathInDir TRUE" }
    relativeInternalPathBothInDir()
        .view{ assert it.text == "relativeInternalPathBothInDirText"; "relativeInternalPathBothInDir TRUE" }
    relativeInternalPathBothInDirButMoved()[1]
        .view{ assert it.text == "relativeInternalPathFileBothInDirButMovedText"; "relativeInternalPathFileBothInDirButMoved TRUE" }
    linkInFolder()
        .view { assert it.resolve("b/linkInFolder.txt").text == "linkInFolderText"; "linkInFolder TRUE" }
}
