package runner;

import algo.needlemanWunsch.NWBanded;
import algo.needlemanWunsch.NWBasic;
import algo.needlemanWunsch.NWGapPenalty;
import algo.sequenceAlignment.ISequenceAlignment;
import algo.smithWaterman.SWBanded;
import algo.smithWaterman.SWBasic;
import org.apache.commons.cli.*;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import util.Common;

import java.io.*;
import java.nio.charset.Charset;
import java.util.*;

import static util.validators.ArgValidator.*;

/**
 * Created by anthony on 08.01.17.
 */
public class Runner {

    private static final int ALGO =    0;
    private static final int INPUT =   1;
    private static final int COUNT =   2;
    private static final int SEQ1 =    3;
    private static final int SEQ2 =    4;
    private static final int SIMILARITY = 5;
    private static final int FEATURE = 6;

    static ImmutablePair<char[], int[][]> aminoPair;

    static Boolean isValid;
    static CommandLine cmd;
    static Options options;
    static List<String> seqList;
    static int count;
    static ISequenceAlignment algo;

    static int penalty;
    static int bandwith;
    static String A;
    static String B;
    static char[] aminoBase;
    static int[][] simMatrix;

    //-a sw-banded -f test.txt -p -5  -sim pam250

    public static void main(String[] args) {

        options = new Options();

        Option optAlgo = new Option(OPT_ALGO_SHORT, OPT_ALGO, true, "algorithm to be used");
        optAlgo.setRequired(true);
        options.addOption(optAlgo);

        Option inputFile = new Option(OPT_FILE_SHORT, OPT_FILE, true, "read input from file [filename]");
        options.addOption(inputFile);

        Option seq1 = new Option(OPT_SEQ1_SHORT, OPT_SEQ1, true, "input sequence 1");
        options.addOption(seq1);

        Option seq2 = new Option(OPT_SEQ2_SHORT, OPT_SEQ2, true, "input sequence 2");
        options.addOption(seq2);

        Option optSim = new Option(OPT_SIM_SHORT, OPT_SIM, true, "similarity matrix to be used");
        optSim.setRequired(true);
        options.addOption(optSim);

        Option pen = new Option(OPT_PEN_SHORT, OPT_PEN, true, "penalty score to be used");
        pen.setRequired(true);
        options.addOption(pen);


        CommandLineParser parser = new DefaultParser();
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {

            System.out.println("Error while parsing args!");
            showHelp();
            System.exit(1);

        }

        seqList = new ArrayList<>();
        count = 2;

        if (cmd.hasOption(OPT_FILE_SHORT)) {

            String fileName = cmd.getOptionValue(OPT_FILE_SHORT);
            seqList = readInputFromFile(fileName);

        } else if (!cmd.hasOption(OPT_SEQ1_SHORT) ||
                   !cmd.hasOption(OPT_SEQ2_SHORT)) {

            System.out.println("Arg '" + OPT_SEQ1 + "' or '" + OPT_SEQ2 + "' error.");
            showHelp();
            System.exit(1);

        } else {

            A = cmd.getOptionValue(OPT_SEQ1_SHORT);
            B = cmd.getOptionValue(OPT_SEQ2_SHORT);

            seqList.add(A);
            seqList.add(B);
        }

        Collections.sort(seqList, seqComparator);

        if (!CollectionUtils.isEmpty(seqList)) {

            penalty = validatePenalty(cmd);
            aminoPair = validateSimilarityMatrix(cmd);

            aminoBase = aminoPair.getLeft();
            simMatrix = aminoPair.getRight();
            bandwith = 3;

            if (!isSeqListValid(seqList, aminoBase)) {
                System.out.println("Some of the provided sequences is invalid against the provided similarity matrix.");
                showHelp();
                System.exit(1);
            }

            A = seqList.get(0);
            B = seqList.get(1);

            algo = validateAlgo(cmd);

            if (null != algo) {

                run();
            }
        }
    }

    private static void run() {

        StringBuilder sb = new StringBuilder();
        sb.append("Computing with provided parameters:\n");
        sb.append("  Algorithm to be used used: \t'").append(getFullName(cmd.getOptionValue(OPT_ALGO_SHORT))).append("'\n");
        sb.append("  For input sequences from ");
        if (cmd.hasOption(OPT_FILE_SHORT)) {
            sb.append("file '").append(cmd.getOptionValue(OPT_FILE_SHORT)).append("'\n");
        } else {
            sb.append("command line\n");
        }
        sb.append("  With penalty score of \t\t'").append(cmd.getOptionValue(OPT_PEN_SHORT)).append("'\n");
        sb.append("  And similarity matrix \t\t'").append(getFullName(cmd.getOptionValue(OPT_SIM_SHORT))).append("'\n");
        System.out.println(sb.toString());

        System.out.println("//==================================================//");
        System.out.println("//===----         Computing Alignment        ----===//");
        System.out.println("//==================================================//");
        System.out.println();

        List<String> res = algo.computeAlignment();

        System.out.println("------ input ------");
        System.out.println(A.toUpperCase());
        System.out.println(B.toUpperCase());

        System.out.println("------ output ------");
        res.forEach(s -> System.out.println(s.toUpperCase()));
        System.out.println();
    }

    private static List<String> readInputFromFile(String fileName) {
        List<String> seqs = new ArrayList<>();
        String line;
        int i = 0;

        try (
                InputStream fis = new FileInputStream(fileName);
                InputStreamReader isr = new InputStreamReader(fis, Charset.forName("UTF-8"));
                BufferedReader br = new BufferedReader(isr);
        ) {
            while (i < count && (line = br.readLine()) != null) {

                seqs.add(line);
            }

        } catch (IOException io) {
            System.out.println("Error reading from file!");
            showHelp();
            System.exit(1);
        }
        return seqs;
    }

    public static ISequenceAlignment validateAlgo(CommandLine cmd) {

        if (cmd.hasOption(OPT_ALGO_SHORT)) {

            String algoStr = cmd.getOptionValue(OPT_ALGO_SHORT).toLowerCase();
            switch (algoStr) {
                case ALGO_NW:
                    return new NWBasic(A, B, penalty, simMatrix, aminoBase);
                case ALGO_NW_BANDED:
                    return new NWBanded(A, B, penalty, bandwith, simMatrix, aminoBase);
                case ALGO_NW_AFFINE:
                    return new NWGapPenalty(A, B, penalty, simMatrix, aminoBase);
                case ALGO_SW:
                    return new SWBasic(A, B, penalty, simMatrix, aminoBase);
                case ALGO_SW_BANDED:
                    return new SWBanded(A, B, penalty, bandwith, simMatrix, aminoBase);

                    default:
                    System.out.println("Arg '" + OPT_ALGO + "' error.");
                    System.out.println("Possible values:");
                    ALGO_VALUES.forEach(System.out::print);
                    System.exit(1);
            }
        } else {
            System.out.println("Arg '" + OPT_ALGO + "' error.");
            showHelp();
            System.exit(1);
        }

        return null;
    }

    public static void showHelp() {
        StringBuilder sb = new StringBuilder();

        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( "biology", options );

        System.out.println();
        sb.append("Possible '" + OPT_ALGO + "' values: \t\t");
        sb.append("[ ");
        ALGO_VALUES.forEach(v -> sb.append(v).append(", "));
        sb.delete(sb.length() - 2, sb.length()).append(" ]");
        System.out.println(sb.toString());

        sb.setLength(0);

        sb.append("Possible '" + OPT_SIM + "' values: \t");
        sb.append("[ ");
        SIM_VALUES.forEach(v -> sb.append(v).append(", "));
        sb.delete(sb.length() - 2, sb.length()).append(" ]");
        System.out.println(sb.toString());

        sb.setLength(0);

        sb.append("Possible '" + OPT_PEN + "' values: \t\t");
        sb.append("< any integer number >");
        System.out.println(sb.toString());

    }

    static Comparator<String> seqComparator = (o1, o2) -> {
        if(o1.length() > o2.length())
            return -1;

        if(o2.length() > o1.length())
            return 1;

        return 0;
    };
}
