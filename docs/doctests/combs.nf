    shapes = Channel.from('circle','square', 'triangle' )

    process combine {
      debug true
      input:
      val shape from shapes
      each color from 'red','blue'
      each size from 1,2,3

      "echo draw $shape $color with size: $size"

    }
