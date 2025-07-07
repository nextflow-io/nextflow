package nextflow.serde.gson

import com.google.gson.TypeAdapter
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter
import nextflow.file.FileHelper

import java.nio.file.Path

class PathAdapter extends TypeAdapter<Path> {
    @Override
    void write(JsonWriter writer, Path value) throws IOException {
        writer.value(value?.toUriString())
    }

    @Override
    Path read(JsonReader reader) throws IOException {
        if (reader.peek() == JsonToken.NULL) {
            reader.nextNull()
            return null
        }
        return FileHelper.asPath(reader.nextString())
    }
}
