    process returnOk {
	echo true
        validExitStatus 0,1,2
         
         script: 
         """
         echo Hello 
         exit 1
	 """    
    }  

