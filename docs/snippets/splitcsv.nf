Channel.of( '10,20,30\n70,80,90' )
    .splitCsv()
    .view { row -> "${row[0]} - ${row[1]} - ${row[2]}" }