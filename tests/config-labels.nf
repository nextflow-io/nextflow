echo true

process alpha {
    /
    echo alpha memry: ${task.memory}
    echo alpha queue: ${task.queue}
    /
}

process beta {
    label 'small'

    /
    echo beta memry: ${task.memory}
    echo beta queue: ${task.queue}
    /
}

process delta {
    label 'big'

    /
    echo delta memry: ${task.memory}
    echo delta queue: ${task.queue}
    /
}

process gamma {
    label 'big'
    memory 4.GB
    queue 'foo'

    /
    echo gamma memry: ${task.memory}
    echo gamma queue: ${task.queue}
    /
}
