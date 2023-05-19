Channel.of( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv( header: ['col1', 'col2', 'col3'], skip: 1 )
    .view { row -> "${row.col1} - ${row.col2} - ${row.col3}" }