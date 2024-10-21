// Example with a JSON object
Channel.of('{"A": 1, "B": [1, 2, 3], "C": {"D": null}}')
    .splitJson()
    .view { v -> "Item: ${v}" }