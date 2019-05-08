package algo.needlemanWunsch;

import java.util.ArrayList;
import java.util.List;

import static util.Common.printCell;
import static util.Common.similarityIndex;
import static util.Const.starFlag;

/**
 * Created by anthony on 15.11.16.
 */
abstract class NWAbstract implements INeedlemanWunsch {

    int[][] m;
    int[][] similarity;
    char[] aminoBases;
    int width, height;
    int penalty;
    String A, B;
    //star
    int score;

    NWAbstract(String A, String B, int penalty, int[][] similarityMatrix, char[] aminoBases) {
        this.A = A;
        this.B = B;
        this.penalty = penalty;
        this.similarity = similarityMatrix;
        this.aminoBases = aminoBases;
    }

    public void setMatrix(int height, int width) {
        this.height = height;
        this.width = width;
        m = new int[height][width];
    }

    abstract void computeMatrix();

    public List<String> computeAlignment() {
        this.score = 0;
        String resA = "";
        String resB = "";

        int i = height - 1;
        int j = width - 1;

        while (i > 0 && j > 0) {

            int score = m[i][j];
            int scoreDiag = m[i - 1][j - 1];
            int scoreUp = m[i][j - 1];
            int scoreLeft = m[i - 1][j];

            //this part is for star alignment
            int idx = similarityIndex(A.charAt(i - 1), aminoBases);
            int jdx = similarityIndex(B.charAt(j - 1), aminoBases);

            int s = (starFlag == idx || starFlag == jdx) ?
                    Integer.MIN_VALUE/2 : this.similarity[idx][jdx];

           /*System.out.println("i/j: " + i + ":" + j + " scr: " + score +
                " scrDiag: " + scoreDiag +
                " scrLeft: " + scoreLeft +
                " s: " + s);*/
            //System.out.println("   A: " + A.charAt(i-1) + " B: " + B.charAt(j-1));


            if (i > 0 && j > 0 && score == scoreDiag + s) {
                this.score += s;
                //System.out.println("a");
                resA = A.charAt(i-1) + resA;
                resB = B.charAt(j-1) + resB;
                i--;
                j--;
            }
            else if (i > 0 && score == scoreLeft + penalty) {
                this.score += penalty;

                resA = A.charAt(i-1) + resA;
                resB = "-" + resB;
                i--;
            }
            else {
                this.score += penalty;

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
