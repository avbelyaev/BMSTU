import java.util.*;

public class Compiler {
    private SortedMap<Position, Message> messages;
    private HashMap<String, Integer> nameCodes;
    private List<String> names;

    public Compiler() {
        messages = new TreeMap<>();
        nameCodes = new HashMap<>();
        names = new ArrayList<>();
    }

    public int addName(String name) {
        if (nameCodes.containsKey(name)) {
            return nameCodes.get(name);
        } else {
            int code = names.size();
            names.add(name);
            nameCodes.put(name, code);
            return code;
        }
    }

    public void addMessage(boolean isErr, Position c, String text) {
        messages.put(c, new Message(isErr, text));
    }

    public void outPutMessages() {
        System.out.println();
        System.out.println("Messages:");

        for (Map.Entry<Position, Message> entry : messages.entrySet()) {
            System.out.print(entry.getValue().isError ? "Error" : "Warning");
            System.out.print(" " + entry.getKey() + ": ");
            System.out.println(entry.getValue().text);
        }
    }

    public MyScanner getScanner(String program) {
        return new MyScanner(program, this);
    }
}
