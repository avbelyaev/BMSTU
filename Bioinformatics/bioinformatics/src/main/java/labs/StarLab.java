package labs;

import algo.multiSequence.StarAlignment;
import algo.sequenceAlignment.ISequenceAlignment;
import util.Const;
import util.Const.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by anthony on 15.02.17.
 */
public class StarLab {

    public static void main(String[] args) {

        String A ="TTGACACCCTCCCAATT";
        String B = "ACCCCAGGCTTTACACAT";
        int penalty = -5;
        int bandWidth = 2;

        /*
S1 YFPHF-DLS-----HGSAQVKAHGKKVG-----DALTLAVAHLDDLPGAL
S2 YFPHF-DLS-----HG-AQVKG—GKKVA-----DALTNAVAHVDDMPNAL
S3 FFPKFKGLTTADQLKKSADVRWHAERII-----NAVNDAVASMDDTEKMS
S4 LFSFLKGTSEVP--QNNPELQAHAGKVFKLVYEAAIQLQVTGVVVTDATL
CO YFPHFKDLS-----HGSAQVKAHGKKVG-----DALTLAVAHVD
         */

        /*
        S1 ACG-TT-GA
S2 ATC-GTCGA
S3 ACGCGA-CC
S4 ACGCGT-TA
         */

        List<String> seqList = new ArrayList<>(Arrays.asList(
                "GATTCA",
                "GTCTGA",
                "GATATT",
                "GTCAGC"
        ));

        ISequenceAlignment sa = new StarAlignment(
                seqList, penalty, Const.BLOSUM45, Const.blosumAminoBases);

        List<String> r = sa.computeAlignment();

        System.out.println("in:");
        System.out.println(A + "\n" + B);
        System.out.println("out:");
        for (String s : r) {
            System.out.println(s);
        }
    }
}
