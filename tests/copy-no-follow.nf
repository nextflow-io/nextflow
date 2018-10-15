process test {
	publishDir "outputDir", mode: 'copyNoFollow'

	output:
	file '*' into testOutput

	"""
	echo "TEST" > testFile.txt
	ln -s testFile.txt testFileLink.txt
	"""
}