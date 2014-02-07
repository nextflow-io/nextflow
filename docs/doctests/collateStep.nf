Channel
  .from(1,2,3,4)
  .collate( 3, 1 )
  .subscribe { println it }
