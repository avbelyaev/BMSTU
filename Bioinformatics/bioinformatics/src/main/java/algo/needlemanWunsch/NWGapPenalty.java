package algo.needlemanWunsch;

import algo.sequenceAlignment.ISequenceAlignment;
import util.Common;
import util.Common.*;

import java.util.ArrayList;
import java.util.List;

import static util.Common.max;
import static util.Common.printCell;
import static util.Common.similarityIndex;

/**
 * Created by anthony on 08.01.17.
 */
public class NWGapPenalty implements INeedlemanWunsch {

    private final char[] aminoBases;

    private int D;  //penalty open
    private int E;  //penalty continuation

    int[][] m;
    int[][] similarityMatrix;
    int width, height;
    String A, B;

    int[][] dist, hor, vert;

    private void setMatrix(int height, int width) {
        this.height = height;
        this.width = width;
        m = new int[height][width];
    }

    public NWGapPenalty(String A, String B, int penalty, int[][] s, char[] aminoBases) {
        this.A = A;
        this.B = B;
        this.similarityMatrix = s;
        this.aminoBases = aminoBases;

        int height = A.length() + 1;
        int width = B.length() + 1;
        this.height = height;
        this.width = width;

        this.D = 12;
        this.E = 9;

        setMatrix(height, width);

        computeMatrix();
    }

    private void computeMatrix() {
        int i, j;

        dist = new int[height + 1][width + 1];
        hor = new int[height + 1][width + 1];
        vert = new int[height + 1][width + 1];

        dist[1][1] = this.similarityMatrix
                [similarityIndex(A.charAt(0), aminoBases)]
                [similarityIndex(B.charAt(0), aminoBases)];

        for (i = 2; i < height; i++) dist[i][1] = this.similarityMatrix
                        [similarityIndex(A.charAt(i-1), aminoBases)]
                        [similarityIndex(B.charAt(0), aminoBases)] - D - E*(i - 2);

        for (i = 2; i < width; i++) dist[1][i] = this.similarityMatrix
                        [similarityIndex(A.charAt(0), aminoBases)]
                        [similarityIndex(B.charAt(i-1), aminoBases)] - D - E*(i - 2);

        for (i = 0; i < height; i++) hor[i][0] = -D - E*(i - 1);
        for (i = 0; i < width; i++) hor[1][i] = -2 * D - E*(i - 1);

        for (i = 0; i < height; i++) vert[i][1] = -2 * D - E*(i - 1);
        for (i = 0; i < width; i++) vert[i][1] = -D - E*(i - 1);


        for (i = 1; i < height - 1; i++) {
            for (j = 1; j < width - 1; j++) {
                dist[i + 1][j + 1] = this.similarityMatrix
                        [similarityIndex(A.charAt(i), aminoBases)]
                        [similarityIndex(B.charAt(j), aminoBases)] + Common.max(dist[i][j], hor[i][j], vert[i][j]);
                hor[i + 1][j] = max(dist[i][j] - D, hor[i][j] - E, vert[i][j] - D);
                vert[i][j + 1] = max(dist[i][j] - D, hor[i][j] - D, vert[i][j] - E);
            }
        }

        for (i = 1; i < height - 1; i++)
            hor[i + 1][width - 1] = max(dist[i][width-1] - D, hor[i][width-1] - E, vert[i][width-1] - D);

        for (i = 1; i < width - 1; i++)
            vert[height-1][i+1] = max(dist[height-1][i] - D, hor[height-1][i] - D, vert[height-1][i] - E);
    }

    public List<String> computeAlignment() {

        String resA = "";
        String resB = "";

        int i = height - 1;
        int j = width - 1;

        while (i > 0 && j > 0) {

            int max = max(dist[i][j], hor[i][j], vert[i][j]);

            if (i > 0 && j > 0 && max == dist[i][j]) {

                //System.out.println("a");
                resA = A.charAt(i-1) + resA;
                resB = B.charAt(j-1) + resB;
                i--;
                j--;
            }
            else if (i > 0 && max == hor[i][j]) {

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
        StringBuilder distStr = new StringBuilder();
        StringBuilder horStr = new StringBuilder();
        StringBuilder vertStr = new StringBuilder();

        int i, j;
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                distStr.append(printCell(dist[i][j])).append(" ");
                horStr.append(printCell(hor[i][j])).append(" ");
                vertStr.append(printCell(vert[i][j])).append(" ");
            }
            distStr.append("\n");
            horStr.append("\n");
            vertStr.append("\n");
        }

        return String.valueOf(distStr) + "\n" +
                horStr + "\n" +
                vertStr + "\n";
    }
}
