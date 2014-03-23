YOURKIT_PORT=${YOURKIT_PORT:-10001}
JVM_OPTIONS="-agentpath:$YOURKIT_HOME/bin/linux-x86-64/libyjpagent.so=port=$YOURKIT_PORT"
exec java $JVM_OPTIONS -jar nextflow "$@"
exit 1
