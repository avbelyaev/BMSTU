import java.io.File;
import java.util.Scanner;

public class complab4 {

    //var 2
    public static void main(String args[]) {

        String text = "";
        Scanner sc;

        try {
            sc = new Scanner(new File("in.txt"));
        } catch (java.io.FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }

        System.out.println("   1234567890123456789");
        int i = 1;
        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            System.out.println(i + " [" + line + "]");
            text += line + "\n";
            i++;
        }

        Compiler compiler = new Compiler();
        MyScanner scanner = new MyScanner(text, compiler);

        System.out.println();
        System.out.println("Tokens:");

        Token t = scanner.nextToken();
        while (null != t) {
            System.out.println(t.toString());
            t = scanner.nextToken();
        }

        scanner.outputComments();

        compiler.outPutMessages();
    }

}
