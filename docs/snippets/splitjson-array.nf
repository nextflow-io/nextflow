// Example with a JSON array
channel.of('[1, null, ["A", {}], true]')
    .splitJson()
    .view { v -> "Item: ${v}" }