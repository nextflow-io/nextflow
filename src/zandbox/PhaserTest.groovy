import java.util.concurrent.Phaser

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

def phaser = new Phaser()

println phaser.register()

println phaser.arrive()
println phaser.arrive()