Channel.of( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv( header: true )
    .view { row -> "${row.alpha} - ${row.beta} - ${row.gamma}" }