Channel.of('[1,null,["A",{}],true]')
    .splitJson()
    .view{"Item: ${it}"}