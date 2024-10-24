Channel.of('{"A": 1, "B": [2, 3, {"C": {"D": null, "E": 4, "F": 5}}]}')
    .splitJson(path: 'B[2].C')
    .view { v -> "Item: ${v}" }