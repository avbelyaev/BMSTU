package algo.needlemanWunsch;

import algo.sequenceAlignment.ISequenceAlignment;
import util.Common;

import java.util.ArrayList;
import java.util.List;

import static util.Common.max;
import static util.Common.printCell;
import static util.Common.similarityIndex;

/**
 * Created by anthony on 15.11.16.
 */
public class NWBanded implements INeedlemanWunsch {

    private final char[] aminoBases;
    private int lenB;
    private int bandWidth;  //k
    private int lenA;        //n

    int[][] m;
    int[][] similarityMatrix;
    int width, height;
    int penalty;
    String A, B;

    void pre(String A, String B, int penalty, int[][] similarityMatrix) {
        this.A = A;
        this.B = B;
        this.penalty = penalty;
        this.similarityMatrix = similarityMatrix;
    }

    public void setMatrix(int height, int width) {
        this.height = height;
        this.width = width;
        m = new int[height][width];
    }

    public NWBanded(String A, String B, int penalty, int bandWidth, int[][] s, char[] aminoBases) {
        pre(A, B, penalty, s);
        this.bandWidth = bandWidth;
        this.lenA = A.length() +1;
        this.lenB = B.length() +1;
        this.aminoBases = aminoBases;

        setMatrix(A.length() + 1, B.length() + 1);

        computeMatrix();
    }

    private boolean insideBand(int i, int j, int k) {
        return -k <= (i - j) && (i - j) <= k;
    }


    void computeMatrix() {
        int i, j, s, h;

        for (i = 0; i < bandWidth; i++)
            m[i][0] = penalty * i;

        for (i = 1; i < bandWidth; i++)
            m[0][i] = penalty * i;

        for (i = 1; i < lenA; i++) {
            for (h = -bandWidth; h <= bandWidth; h++) {

                j = i + h;

                if (1 <= j && j < lenB) {

                    //System.out.println("i/j: " + i + ":" + j);
                    s = this.similarityMatrix
                            [Common.similarityIndex(A.charAt(i-1), aminoBases)]
                            [Common.similarityIndex(B.charAt(j-1), aminoBases)];
                    m[i][j] = m[i - 1][j - 1] + s;

                    if (insideBand(i - 1, j, bandWidth)) {

                        int delete = m[i - 1][j] + penalty;
                        m[i][j] = max( m[i][j], delete);
                    }

                    if (insideBand(i, j - 1, bandWidth)) {

                        int insert = m[i][j - 1] + penalty;
                        m[i][j] = max( m[i][j], insert);
                    }
                }
            }
        }
    }

    public List<String> computeAlignment() {
        int s;

        String resA = "";
        String resB = "";

        int i = height - 1;
        int j = width - 1;

        while (i > 0 && j > 0) {

            int score = m[i][j];
            int scoreDiag = m[i - 1][j - 1];
            int scoreUp = m[i][j - 1];
            int scoreLeft = m[i - 1][j];

            s = this.similarityMatrix
                    [Common.similarityIndex(A.charAt(i-1), aminoBases)]
                    [Common.similarityIndex(B.charAt(j-1), aminoBases)];
            /*System.out.println("i/j: " + i + ":" + j + " scr: " + score +
                    " scrDiag: " + scoreDiag +
                    " scrLeft: " + scoreLeft +
                    " s: " + s);
            System.out.println("   A: " + A.charAt(i-1) + " B: " + B.charAt(j-1));
            System.out.println("   res A:" + resA + " B: " + resB);
            */

            if (i > 0 && j > 0 && score == scoreDiag + s) {

                //System.out.println("diag");
                resA = A.charAt(i-1) + resA;
                resB = B.charAt(j-1) + resB;
                i--;
                j--;
            }
            else if (i > 0 && score == scoreLeft + penalty) {

                //System.out.println("left");
                resA = A.charAt(i-1) + resA;
                resB = "-" + resB;
                i--;
            }
            else {

                //System.out.println("top");
                resA = "-" + resA;
                resB = B.charAt(j-1) + resB;
                j--;
            }
        }

        while (i > 0) {
            resA = A.charAt(i-1) + resA;
            resB = "-" + resB;
            i--;
        }

        while (j > 0) {
            resA = "-" + resA;
            resB = B.charAt(j-1) + resB;
            j--;
        }

        List<String> res = new ArrayList<String>();
        res.add(resA);
        res.add(resB);

        return res;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();

        int i, j;
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                s.append(printCell(m[i][j])).append(" ");
            }
            s.append("\n");
        }

        return s.toString();
    }
}
