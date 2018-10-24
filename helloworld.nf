process sayHello {
	
	output:
   	file 'result.txt' into numbers

   	'''
   	echo $RANDOM > result.txt
   	'''	
       
}

process sayGoodbye {
	"""
	printf 'Goodbye!'
	"""
}

