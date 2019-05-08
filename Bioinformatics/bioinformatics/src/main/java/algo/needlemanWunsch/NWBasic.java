package algo.needlemanWunsch;

import util.Const;

import static util.Common.max;
import static util.Common.similarityIndex;
import static util.Const.starFlag;

/**
 * Created by anthony on 27.10.16.
 */
public class NWBasic extends NWAbstract implements INeedlemanWunsch {

    public NWBasic(String A, String B, int penalty, int[][] s, char[] aminoBases) {
        super(A, B, penalty, s, aminoBases);

        int height = A.length() + 1;
        int width = B.length() + 1;

        setMatrix(height, width);

        computeMatrix();
    }

    @Override
    void computeMatrix() {
        int i, j;

        for (i = 0; i < height; i++)
            m[i][0] = penalty * i;

        for (i = 1; i < width; i++)
            m[0][i] = penalty * i;

        for (i = 1; i < height; i++) {
            for (j = 1; j < width; j++) {

             //System.out.println("a/b: " + A.charAt(i - 1) + "/" + B.charAt(j - 1));

                //this part is for star alignment
                int idx = similarityIndex(A.charAt(i - 1), aminoBases);
                int jdx = similarityIndex(B.charAt(j - 1), aminoBases);

                int s = (starFlag == idx || starFlag == jdx) ?
                        Integer.MIN_VALUE /2 : this.similarity[idx][jdx];

                int match = m[i - 1][j - 1] + s;
                int delete = m[i - 1][j] + penalty;
                int insert = m[i][j - 1] + penalty;

                m[i][j] = max(match, insert, delete);
            }
        }
    }

    public int getScore() {
        return score;
    }
}
