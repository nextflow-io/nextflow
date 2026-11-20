import java.io.*;
import java.lang.classfile.*;
import java.lang.classfile.constantpool.*;
import java.lang.reflect.AccessFlag;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.*;
import java.util.*;

/**
 * Reports linkage errors of compiled classes against a classpath: referenced classes, methods and fields
 * that do not exist, and inherited abstract methods with no implementation.
 * Usage: java Linkage.java <classpath> <classes-dir>...
 * Exits with 3 when it finds problems, so other failures can be told apart.
 */
public class Linkage {
    static final Map<String, Path> local = new HashMap<>();
    static final Map<String, Optional<ClassModel>> cache = new HashMap<>();
    // Reads class bytes only: jars on -cp would make the source launcher load their services.
    static ClassLoader classpath;

    public static void main(String[] args) throws IOException {
        var urls = new ArrayList<URL>();
        for (String entry : args[0].split(File.pathSeparator))
            if (!entry.isEmpty()) urls.add(Path.of(entry).toUri().toURL());
        classpath = new URLClassLoader(urls.toArray(URL[]::new), ClassLoader.getPlatformClassLoader());
        for (String arg : Arrays.copyOfRange(args, 1, args.length)) {
            Path root = Path.of(arg);
            try (var files = Files.walk(root)) {
                files.filter(p -> p.toString().endsWith(".class")).forEach(p -> {
                    String name = root.relativize(p).toString().replace(File.separatorChar, '/');
                    local.put(name.substring(0, name.length() - 6), p);
                });
            }
        }
        var problems = new TreeSet<String>();
        for (String name : local.keySet()) {
            ClassModel cm = load(name);
            for (PoolEntry e : cm.constantPool()) {
                String problem = switch (e) {
                    case MemberRefEntry m -> checkMember(m);
                    case ClassEntry c -> checked(c.asInternalName()) && load(c.asInternalName()) == null ? "missing class " + c.asInternalName() : null;
                    default -> null;
                };
                if (problem != null) problems.add(name + ": " + problem);
            }
            if (!cm.flags().has(AccessFlag.ABSTRACT) && !cm.flags().has(AccessFlag.INTERFACE))
                for (String m : unimplemented(cm)) problems.add(name + ": abstract method not implemented " + m);
        }
        problems.forEach(System.out::println);
        System.exit(problems.isEmpty() ? 0 : 3);
    }

    /** Leaves out array and primitive types. */
    static boolean checked(String name) {
        return !name.startsWith("[") && name.contains("/");
    }

    static ClassModel load(String name) {
        return cache.computeIfAbsent(name, n -> {
            try {
                Path p = local.get(n);
                if (p != null) return Optional.of(ClassFile.of().parse(p));
                try (InputStream in = classpath.getResourceAsStream(n + ".class")) {
                    return in == null ? Optional.empty() : Optional.of(ClassFile.of().parse(in.readAllBytes()));
                }
            } catch (IOException e) { throw new UncheckedIOException(e); }
        }).orElse(null);
    }

    /** All supertypes of a class, itself included; null when one of them is not on the classpath. */
    static List<ClassModel> supertypes(ClassModel cm) {
        var result = new ArrayList<ClassModel>();
        var seen = new HashSet<String>();
        var todo = new ArrayDeque<ClassModel>(List.of(cm));
        while (!todo.isEmpty()) {
            ClassModel c = todo.pop();
            if (!seen.add(c.thisClass().asInternalName())) continue;
            result.add(c);
            var parents = new ArrayList<>(c.interfaces());
            c.superclass().ifPresent(parents::add);
            for (ClassEntry parent : parents) {
                ClassModel p = load(parent.asInternalName());
                if (p == null) return null;
                todo.add(p);
            }
        }
        return result;
    }

    static String checkMember(MemberRefEntry m) {
        String owner = m.owner().asInternalName();
        if (!checked(owner)) return null;
        ClassModel cm = load(owner);
        if (cm == null) return null;  // reported as a missing class
        List<ClassModel> types = supertypes(cm);
        if (types == null) return null;  // can't tell when a supertype is a dependency not on the classpath
        String name = m.name().stringValue(), desc = m.type().stringValue();
        boolean found = types.stream().anyMatch(t -> m instanceof FieldRefEntry
            ? t.fields().stream().anyMatch(f -> f.fieldName().equalsString(name) && f.fieldType().equalsString(desc))
            : t.methods().stream().anyMatch(x -> x.methodName().equalsString(name) && x.methodType().equalsString(desc)));
        return found ? null : "missing " + (m instanceof FieldRefEntry ? "field " : "method ") + owner + "." + name + desc;
    }

    /** A method that can implement or override: private and static ones cannot. */
    static MethodModel find(ClassModel c, String sig) {
        return c.methods().stream()
            .filter(m -> !m.flags().has(AccessFlag.STATIC) && !m.flags().has(AccessFlag.PRIVATE) && sig.equals(m.methodName().stringValue() + m.methodType().stringValue()))
            .findFirst().orElse(null);
    }

    /** Abstract methods that JVM method selection (JVMS 5.4.6) would not find an implementation for. */
    static List<String> unimplemented(ClassModel cm) {
        List<ClassModel> types = supertypes(cm);
        if (types == null) return List.of();
        var abstracts = new TreeMap<String, String>();
        for (ClassModel t : types)
            for (MethodModel m : t.methods())
                if (m.flags().has(AccessFlag.ABSTRACT))
                    abstracts.putIfAbsent(m.methodName().stringValue() + m.methodType().stringValue(), t.thisClass().asInternalName());
        var result = new ArrayList<String>();
        for (var e : abstracts.entrySet())
            if (!implemented(cm, types, e.getKey())) result.add(e.getValue() + "." + e.getKey());
        return result;
    }

    static boolean implemented(ClassModel cm, List<ClassModel> types, String sig) {
        // The first declaration in the superclass chain wins.
        for (ClassModel c = cm; c != null; c = c.superclass().map(s -> load(s.asInternalName())).orElse(null)) {
            MethodModel m = find(c, sig);
            if (m != null) return !m.flags().has(AccessFlag.ABSTRACT);
        }
        // Otherwise the maximally-specific interface declaration: a sub-interface that redeclares it abstract hides a default.
        var declaring = types.stream().filter(t -> t.flags().has(AccessFlag.INTERFACE) && find(t, sig) != null).toList();
        var specific = declaring.stream().filter(i -> declaring.stream()
            .noneMatch(j -> j != i && supertypes(j).contains(i))).toList();
        return specific.stream().filter(i -> !find(i, sig).flags().has(AccessFlag.ABSTRACT)).count() == 1;
    }
}
