package util.validators;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.lang3.tuple.ImmutablePair;
import runner.Runner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static util.Const.*;
import static util.Const.DNAFull;
import static util.Const.dnaFullAminoBases;

/**
 * Created by anthony on 18.02.17.
 */
public class ArgValidator {

    //args
    public static final String OPT_ALGO_SHORT = "a";
    public static final String OPT_FILE_SHORT = "f";
    public static final String OPT_SEQ1_SHORT = "s1";
    public static final String OPT_SEQ2_SHORT = "s2";
    public static final String OPT_SIM_SHORT = "sim";
    public static final String OPT_PEN_SHORT = "p";

    public static final String OPT_ALGO = "algo";
    public static final String OPT_FILE = "file";
    public static final String OPT_SEQ1 = "seq1";
    public static final String OPT_SEQ2 = "seq2";
    public static final String OPT_SIM = "sim-matrix";
    public static final String OPT_PEN = "penalty";


    //sim matrices
    public static final String SIM_BASIC = "basic";
    public static final String SIM_SIMPLE = "simple";
    public static final String SIM_BLOSUM45 = "blosum45";
    public static final String SIM_BLOSUM62 = "blosum62";
    public static final String SIM_DNAFULL = "dnafull";
    public static final String SIM_PAM250 = "pam250";

    public static final List<String> SIM_VALUES = new ArrayList<>(Arrays.asList(
            SIM_BASIC,
            SIM_SIMPLE,
            SIM_BLOSUM45,
            SIM_BLOSUM62,
            SIM_DNAFULL,
            SIM_PAM250
    ));

    public static final String SIM_NAME_BASIC = "Basic similarity matrix";
    public static final String SIM_NAME_SIMPLE = "Simple similarity matrix";
    public static final String SIM_NAME_BLOSUM45 = "BLOSUM 45";
    public static final String SIM_NAME_BLOSUM62 = "BLOSUM 62";
    public static final String SIM_NAME_DNAFULL = "DNA Full";
    public static final String SIM_NAME_PAM250 = "PAM 250";


    //algos
    public static final String ALGO_NW = "nw";
    public static final String ALGO_NW_BANDED = "nw-banded";
    public static final String ALGO_NW_AFFINE = "nw-affine";
    public static final String ALGO_SW = "sw";
    public static final String ALGO_SW_BANDED = "sw-banded";

    public static final List<String> ALGO_VALUES = new ArrayList<>(Arrays.asList(
            ALGO_NW,
            ALGO_NW_BANDED,
            ALGO_NW_AFFINE,
            ALGO_SW,
            ALGO_SW_BANDED
    ));

    public static final String ALGO_NAME_NW = "Needleman-Wunsch Basic";
    public static final String ALGO_NAME_NW_BANDED = "Needleman-Wunsch Banded";
    public static final String ALGO_NAME_NW_AFFINE = "Needleman-Wunsch Affine Gap Penalty";
    public static final String ALGO_NAME_SW = "Smith-Waterman Basic";
    public static final String ALGO_NAME_SW_BANDED = "Swmith-Waterman Banded";

    public static ImmutablePair<char[], int[][]> validateSimilarityMatrix(CommandLine cmd) {
        ImmutablePair<char[], int[][]> aminoPair = null;

        if (cmd.hasOption(OPT_SIM_SHORT)) {

            String sim = cmd.getOptionValue(OPT_SIM_SHORT).toLowerCase();

            switch (sim) {
                case SIM_BASIC:
                    return new ImmutablePair<>(basicAmeno, basicSimilarity);
                case SIM_SIMPLE:
                    return new ImmutablePair<>(basicAmeno, simpleMatrix(-5));
                case SIM_BLOSUM45:
                    return new ImmutablePair<>(blosumAminoBases, BLOSUM45);
                case SIM_BLOSUM62:
                    return new ImmutablePair<>(blosumAminoBases, BLOSUM62);
                case SIM_DNAFULL:
                    return new ImmutablePair<>(dnaFullAminoBases, DNAFull);
                case SIM_PAM250:
                    return new ImmutablePair<>(pam250AminoBases, PAM250);
                default:
                    System.out.println("Arg '" + OPT_SIM + "' error.");
                    Runner.showHelp();
                    System.exit(1);
            }
        }
        return aminoPair;
    }

    public static int validatePenalty(CommandLine cmd) {
        int pen = penaltyDefaultValue;
        if (cmd.hasOption(OPT_PEN_SHORT)) {

            String penStr = cmd.getOptionValue(OPT_PEN_SHORT);
            try {
                pen = Integer.valueOf(penStr);

            } catch (Exception e) {
                System.out.println("Arg '" + OPT_PEN + "' error.");
                Runner.showHelp();
                System.exit(1);
            }
        }
        return pen;
    }

    public static boolean isValidSequence(String sequence, char[] charset) {

        String chars = new String(charset);
        chars = "[" + chars + "]*";

        return sequence.matches(chars);
    }

    public static boolean isSeqListValid(List<String> seqList, char[] charset) {

        for (String s : seqList) {
            if (!isValidSequence(s, charset))
                return false;
        }
        return true;
    }

    public static String getFullName(String shortName) {
        switch (shortName) {
            case ALGO_NW:
                return ALGO_NAME_NW;
            case ALGO_NW_BANDED:
                return ALGO_NAME_NW_BANDED;
            case ALGO_NW_AFFINE:
                return ALGO_NAME_NW_AFFINE;
            case ALGO_SW:
                return ALGO_NAME_SW;
            case ALGO_SW_BANDED:
                return ALGO_NAME_SW_BANDED;

            case SIM_BASIC:
                return SIM_NAME_BASIC;
            case SIM_SIMPLE:
                return SIM_NAME_SIMPLE;
            case SIM_BLOSUM45:
                return SIM_NAME_BLOSUM45;
            case SIM_BLOSUM62:
                return SIM_NAME_BLOSUM62;
            case SIM_DNAFULL:
                return SIM_NAME_DNAFULL;
            case SIM_PAM250:
                return SIM_NAME_PAM250;

            default:
                return "Unknown'. App error possible!";
        }
    }
}
