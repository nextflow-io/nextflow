     tuple = Channel.from( [1, 'alpha', 'Lorem ipsum dolor sit amet'], [2, 'beta', 'consectetur adipiscing elit'] )

     process setExample {
         input:
         set val(x), env(VAR_X), file('latin.txt')  from tuple 

         """
         echo Processing $x
         cat latin.txt
         """

     }

