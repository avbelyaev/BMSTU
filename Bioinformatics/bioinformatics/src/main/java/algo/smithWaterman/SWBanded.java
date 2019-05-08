package algo.smithWaterman;

import util.Common;

import java.util.ArrayList;
import java.util.List;

import static util.Common.max;
import static util.Common.similarityIndex;

/**
 * Created by anthony on 05.02.17.
 */
public class SWBanded implements ISmithWaterman {

    private final char[] aminoBases;
    private int lenB;
    private int bandWidth;  //k
    private int lenA;        //n

    int[][] m;
    int[][] similarityMatrix;
    int width, height;
    int penalty;
    String A, B;

    int max, maxi, maxj;

    private void setMatrix(int height, int width) {
        this.height = height;
        this.width = width;
        m = new int[height][width];
    }

    public SWBanded(String A, String B, int penalty, int bandwidth, int[][] s, char[] aminoBases) {
        this.A = A;
        this.B = B;
        this.penalty = penalty;
        this.similarityMatrix = s;
        this.aminoBases = aminoBases;

        this.bandWidth = bandwidth;
        int height = A.length() + 1;
        int width = B.length() + 1;

        setMatrix(height, width);

        computeMatrix();
    }

    private boolean insideBand(int i, int j, int k) {
        return -k <= (i - j) && (i - j) <= k;
    }

    private void computeMatrix() {
        int i, j, h, s;

        max = 0;

        for (i = 0; i < bandWidth; i++) {
            m[i][0] = 0;
            m[0][i] = 0;
        }

        for (i = 1; i < height; i++) {
            for (h = -bandWidth; h <= bandWidth; h++) {

                j = i + h;

                if (1 <= j && j < width) {

                    s = this.similarityMatrix
                            [Common.similarityIndex(A.charAt(i-1), aminoBases)]
                            [Common.similarityIndex(B.charAt(j-1), aminoBases)];

                    int match = m[i - 1][j - 1] + s;
                    int delete = m[i - 1][j] + penalty;
                    int insert = m[i][j - 1] + penalty;

                    m[i][j] = max(0, match, insert, delete);

                    if (m[i][j] > max) {
                        max = m[i][j];
                        maxi = i;
                        maxj = j;
                    }

                    if (insideBand(i - 1, j, bandWidth)) {

                        //System.out.println("a");
                        m[i][j] = max(0, m[i][j], delete);
                    }

                    if (insideBand(i, j - 1, bandWidth)) {

                        //System.out.println("b");
                        m[i][j] = max(0, m[i][j], insert);
                    }
                }
            }
        }
    }

    public List<String> computeAlignment() {

        String resA = "";
        String resB = "";

        int j = width - 1;
        int i = height - 1;

        while (i > maxi) {
            resA = A.charAt(i-1) + resA;
            resB = "-" + resB;
            i--;
        }

        while (j > maxj) {
            resA = "-" + resA;
            resB = B.charAt(j-1) + resB;
            j--;
        }

        //System.out.println("max[" + maxi + "][" + maxj + "]: " + max);
        i = maxi;
        j = maxj;

        while (i > 0 && j > 0 && m[i][j] > 0) {

            int score = m[i][j];
            int scoreDiag = m[i - 1][j - 1];
            int scoreUp = m[i][j - 1];
            int scoreLeft = m[i - 1][j];

            int s = this.similarityMatrix
                    [Common.similarityIndex(A.charAt(i-1), aminoBases)]
                    [Common.similarityIndex(B.charAt(j-1), aminoBases)];

            /*System.out.println("i/j: " + i + ":" + j + " scr: " + score +
                    " scrDiag: " + scoreDiag +
                    " scrLeft: " + scoreLeft +
                    " s: " + s);
            System.out.println("   A: " + A.charAt(i-1) + " B: " + B.charAt(j-1));
            */

            if (i > 0 && j > 0 && score == scoreDiag + s) {

                //System.out.println("a");
                resA = A.charAt(i-1) + resA;
                resB = B.charAt(j-1) + resB;
                i--;
                j--;
            }
            else if (i > 0 && score == scoreLeft + penalty) {


                resA = A.charAt(i-1) + resA;
                resB = "-" + resB;
                i--;
            }
            else {

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
                s.append(Common.printCell(m[i][j]));
                if (i == maxi && j == maxj)
                    s.append("X");
                s.append(" ");
            }
            s.append("\n");
        }

        return s.toString();
    }

}
