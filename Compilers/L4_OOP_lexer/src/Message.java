public class Message {
    public final Boolean isError;
    public final String text;

    public Message(Boolean isError, String text) {
        this.isError = isError;
        this.text = text;
    }
}
