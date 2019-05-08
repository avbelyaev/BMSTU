package algo.multiSequence;

import algo.needlemanWunsch.NWBasic;
import algo.sequenceAlignment.ISequenceAlignment;
import org.apache.commons.collections.CollectionUtils;
import util.Const;

import java.util.*;

/**
 * Created by anthony on 08.01.17.
 */
public class StarAlignment implements ISequenceAlignment{

    private final int[][] similarityMatrix;
    private final int penalty;
    private final char[] aminoBases;

    List<String> seqList;
    List<String> alreadyAligned;

    public StarAlignment(
            List<String> seqList, int penalty, int[][] similarityMatrix, char[] aminoBases) {
        this.seqList = seqList;
        this.penalty = penalty;
        this.similarityMatrix = similarityMatrix;
        this.aminoBases = aminoBases;

        this.seqList = seqList;
    }

    public List<String> computeAlignment() {
        int i;
        int n = seqList.size();

        int centerIndex = findCenter(seqList);
        String Scenter = seqList.get(centerIndex);

        alreadyAligned = new LinkedList<>();

        Stack<MyPair> stack = new Stack<>();
        List<String> tmpInit = new ArrayList<>();

        for (i = 1; i < n; i++) {
            if (i != centerIndex) {

                System.out.println("i: " + i);

                //-4-
                String Si = seqList.get(i);

                System.out.println("Si:       " + Si);
                System.out.println("Scenter:  " + Scenter);

                ISequenceAlignment sw = new NWBasic(Si, Scenter, penalty, Const.basicSimilarity, Const.basicAmeno);
                List<String> res = sw.computeAlignment();
                System.out.println("----------------------------");

                //-5-
                String Si_new = res.get(0);
                String Scenter_new = res.get(1);

                System.out.println("Si_new:   " + Si_new);
                System.out.println("Sctr_new: " + Scenter_new);


                tmpInit.add(Si_new);

                stack.push(new MyPair(Scenter_new, new ArrayList<>(tmpInit)));

               /* List<Integer> gapIndices = findGaps(Scenter_new);
                alreadyAligned.add(Si_new);
                adjustGaps(gapIndices, alreadyAligned);

                System.out.println("already done:");
                alreadyAligned.forEach(System.out::println);
                //System.out.println("\n");
                //-7-
                Scenter = Scenter_new;*/
            }


            /*
                "GATTCA",
                "GTCTGA",
                "GATATT",
                "GTCAGC"
             */
        }

        while (stack.size() > 1) {

            MyPair seq1 = stack.peek();
            stack.pop();
            MyPair seq2 = stack.peek();
            stack.pop();

            MyPair merged = merge(seq1, seq2);

            stack.push(merged);
        }

        System.out.println("res:");
        System.out.println(stack.peek().first);
        stack.peek().second.forEach(System.out::println);
        //System.out.println(Scenter);
        /*List<String> retList = new ArrayList<>();
        retList.add(seqList.get(0));
        retList.retainAll(alreadyAligned);
        */

        return new ArrayList<>();
    }

    private MyPair merge(MyPair a, MyPair b) {

        System.out.println("\nmerge{..}");
        System.out.println("a:");
        System.out.println(a.first);
        a.second.forEach(System.out::println);
        System.out.println("b:");
        System.out.println(b.first);
        b.second.forEach(System.out::println);

        int i = 0;
        List<String> resVect = new ArrayList<>();
        StringBuilder res = new StringBuilder();

        while ( i < a.first.length() && i < b.first.length()) {

            if (a.first.charAt(0) == b.first.charAt(0)) {

                res.append(a.first.charAt(i));
                i++;
            } else {

                if (a.first.charAt(i) == '-') {

                    res.append("-");
                    b.first = b.first.substring(0, i) + "-" + b.first.substring(i, b.first.length() - 1);

                    int j = 0;
                    for (String str : b.second) {
                        b.second.set(j, (str.substring(0, i) + "-" + str.substring(i, str.length() - 1)));
                        j++;
                    }
                    i++;

                } else {

                    res.append("-");
                    a.first = a.first.substring(0, i) + "-" + a.first.substring(i, a.first.length() - 1);

                    int j = 0;
                    for (String str : a.second) {
                        a.second.set(j, (str.substring(0, i) + "-" + str.substring(i, str.length() - 1)));
                        j++;
                    }
                    i++;
                }
            }
        }

        for (String str : b.second) {
            resVect.add(str);
        }
        for (String str : a.second) {
            resVect.add(str);
        }

        System.out.println("res:");
        System.out.println(res.toString());
        resVect.forEach(System.out::println);
        System.out.println();

        return new MyPair(res.toString(), resVect);
    }

    private static class MyPair {
        String first;
        List<String> second;

        public MyPair(String first, List<String> second) {
            this.first = first;
            this.second = second;
        }

        public String getFirst() {
            return first;
        }

        public void setFirst(String first) {
            this.first = first;
        }

        public List<String> getSecond() {
            return second;
        }

        public void setSecond(List<String> second) {
            this.second = second;
        }
    }

    private List<Integer> findGaps(String str) {
        List<Integer> gaps = new ArrayList<>();

        for (int i = 0; i < str.length(); i++) {
            if ("-".equals(String.valueOf(str.charAt(i)))) {
                gaps.add(i);
                System.out.println("gap[" + i + "]");
            }
        }
        return gaps;
    }

    private void adjustGaps(final List<Integer> gapIndices, List<String> seqs) {
        System.out.print("adjustGaps{..}");
        if (!CollectionUtils.isEmpty(gapIndices)) {

            int i = 0;
            for (String s : seqs) {
                seqs.set(i++, insertGap(s, gapIndices));
            }

        } else {
            System.out.println(" - skipping");
        }
    }

    private String insertGap(String str, final List<Integer> gapIndices) {
        System.out.println("insertGap{..}");
        List<Integer> gapsToBeAdded = new ArrayList<>();

        System.out.println("str[" + str.length() + "]: " + str);

        /*int j = 0;
        for (Integer gap : gapIndices) {

            if (gap < str.length() && !"-".equals(String.valueOf(str.charAt(gap)))) {
                gapsToBeAdded.add(gap);
                j++;
            }
        }
        while (j < gapIndices.size()) {
            gapsToBeAdded.add(gapIndices.get(j));
            j++;
        }*/
        gapIndices.forEach(gap -> {
            if (gap < str.length() && !"-".equals(String.valueOf(str.charAt(gap)))) {
                gapsToBeAdded.add(gap);
            }
        });

        String strCopy = str;
        StringBuilder s = new StringBuilder();

        gapsToBeAdded.sort(Collections.reverseOrder());

        System.out.println("gaps:");
        gapsToBeAdded.forEach(System.out::println);

        int i = 0;
        for (Integer gap : gapsToBeAdded) {

            String pre = strCopy.substring(0, gap + i);
            String post = strCopy.substring(gap + i, strCopy.length() - 1);
            System.out.println("pre: " + pre + ", post: " + post);
            s.append(pre).append("-").append(post);

            strCopy = s.toString();
            //i++;
        }

        System.out.println("new s: " + s);
        return s.toString();
    }

    private void swapSeq(int i, int j, List<String> seqs) {
        String tmp = seqs.get(i);
        seqs.set(i, seqs.get(j));
        seqs.set(j, tmp);
    }

    private int findCenter(List<String> seqList) {
        int n = seqList.size();
        int[] score = new int[n];
        int i, j;

        for (i = 0; i < n - 1; i++) {
            for (j = i + 1; j < n ;j++) {

                NWBasic nw = new NWBasic(
                        seqList.get(i),
                        seqList.get(j),
                        penalty,
                        Const.basicSimilarity,
                        Const.basicAmeno
                );
                nw.computeAlignment();

                score[i] += nw.getScore();
                score[j] += nw.getScore();
            }
        }

        int maxScore = Integer.MIN_VALUE;
        j = -1;

        for (i = 0; i < n ; i++) {
            if (score[i] > maxScore)
                maxScore = score[j = i];
        }

        return j;
    }
}
