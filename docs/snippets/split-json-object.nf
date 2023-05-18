Channel.of('{"A":1,"B":[1,2,3],"C":{"D":null}}')
    .splitJson()
    .view{"Item: ${it}"}