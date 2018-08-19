//hadlab8: Trident Storm
Лабораторная работа 8. Разработка простой топологии trident.

Задача:
-Требуется разработать топологию для анализа данных перелетов — получить
	распределение количества опоздавших или отмененных рейсов по дням недели.
-а.Требуется разработать Spout который будет опрашивать директорию и в
	случае обнаружения в ней новых файлов читать ее построчно и генерировать
	batch, cостоящий из всех строк файла. После окончания файла требуется
	перенести файл в папку для обработанных файлов
-б. Разрабатываем Function которая будет принимать строку и разбивать ее на
	tuple состоящий из следующих данных 664600583_T_ONTIME_sample.csv :
	день недели, время опоздания, признак совершенного рейса
-в. Разрабатываем фильтр отбрасывающий записи рейсов которые прилетели
	вовремя
-г. Разрабатываем CombinerAggregator который принимает tuple с данными
	совершенного рейса и подсчитывает аггрегированные данные <день
	недели>:<количество рейсов>
-д. Разрабатываем фильтр который будет принимать результат с данными
	CombinerAggregator и печатать их в консоль



//запуск:
//mnv package
//mvn exec:java -Dexec.mainClass="TridentLab"



//TridentLab.java
import backtype.storm.Config;
import backtype.storm.LocalCluster;
import backtype.storm.topology.TopologyBuilder;
import backtype.storm.tuple.Fields;
import storm.trident.TridentTopology;

public class TridentLab {

    public static void main(String[] args) {

        TridentTopology topology = new TridentTopology();

        topology.newStream("generator", new BatchSpout())
                .parallelismHint(1)
                .shuffle()
                .each(new Fields("sentence"), new SplitFunction(), new Fields("day", "delay", "cancellation"))
                .each(new Fields("delay", "cancellation"), new DelayFilter())
                .groupBy(new Fields("day"))
                .aggregate(new Fields("day"), new DayAggregator(), new Fields("count"))
                .parallelismHint(7)
                .each(new Fields("day", "count"), new PrinterFunction(), new Fields());

        Config config = new Config();

        LocalCluster cluster = new LocalCluster();
        cluster.submitTopology("poll", config, topology.build());
    }

}




//BatchSpout.java
import backtype.storm.task.TopologyContext;
import backtype.storm.tuple.Fields;
import backtype.storm.tuple.Values;
import com.google.common.base.Charsets;
import com.google.common.io.Files;
import storm.trident.operation.TridentCollector;
import storm.trident.spout.IBatchSpout;

import java.io.*;
import java.lang.Override;
import java.util.Map;

public class BatchSpout implements IBatchSpout {
    private TopologyContext context;
    private File dir;

    public  BatchSpout() {

//2 directories should exist: in_files, out_files
//in_files must contain 664600583_T_ONTIME_sample.csv

        dir = new File("/home/anthony/hadlab8/in_files");
    }

    public Fields getOutputFields() {
        return new Fields("sentence");
    }

    public void open(Map map, TopologyContext topologyContext) {
        context = topologyContext;
    }

    public void emitBatch(long l, TridentCollector tridentCollector) {

        File files[] = dir.listFiles();

        if (null != files && 0 != files.length) {

            File curr_file = files[0];
            try {
                BufferedReader reader = new BufferedReader(new InputStreamReader(
                        new FileInputStream(curr_file), Charsets.UTF_8));

                while (true) {

                    String line = reader.readLine();
                    if (null == line) {
                        break;
                    }

                    tridentCollector.emit(new Values(line));
                }

                reader.close();
                Files.move(curr_file, new File("/home/anthony/hadlab8/out_files/" + curr_file.getName()));

            } catch (FileNotFoundException e1) { e1.printStackTrace(); }
              catch (IOException e2) { e2.printStackTrace(); }

        }
    }

    public void ack(long l) {}

    public void close() {}

    public Map getComponentConfiguration() { return null; }

}




//SplitFunction.java
import backtype.storm.tuple.Values;
import storm.trident.operation.BaseFunction;
import storm.trident.operation.TridentCollector;
import storm.trident.tuple.TridentTuple;

public class SplitFunction extends BaseFunction {
    private static final int DAY_OF_WEEK = 4;
    private static final int ARR_DELAY_NEW = 18;
    private static final int CANCELLED = 19;

    public void execute(TridentTuple tuple, TridentCollector collector) {
        String[] columns = tuple.getString(0).split(",");
        if (!"\"YEAR\"".equals(columns[0])) {

            collector.emit(new Values(
                    columns[DAY_OF_WEEK],
                    columns[ARR_DELAY_NEW],
                    columns[CANCELLED]));
        }
    }

}




//DelayFilter.java
import storm.trident.operation.BaseFilter;
import storm.trident.tuple.TridentTuple;

public class DelayFilter extends BaseFilter {
    private static final float ON_TIME = (float)0;

    public boolean isKeep(TridentTuple tuple) {
        String delay = tuple.getString(0);
        String cancellation = tuple.getString(1);

        if (null != cancellation && !cancellation.equals("") && cancellation.equals("1.00")) {
            return true;
        } else {
            if (null != delay && !delay.equals("") && ON_TIME <= Float.parseFloat(delay)) {
                return true;
            }
        }
        return false;
    }
}




//DayAggregator.java
import storm.trident.operation.CombinerAggregator;
import storm.trident.tuple.TridentTuple;

public class DayAggregator implements CombinerAggregator<Long> {

    public Long init(TridentTuple tridentTuple) {
        return 1L;
    }

    public Long combine(Long a, Long b) {
        return a + b;
    }

    public Long zero() {
        return 0L;
    }
}




//PrinterFunction.java
import storm.trident.operation.BaseFunction;
import storm.trident.operation.TridentCollector;
import storm.trident.tuple.TridentTuple;

public class PrinterFunction extends BaseFunction {
    public void execute(TridentTuple tridentTuple, TridentCollector tridentCollector) {
        if (!tridentTuple.toString().equals("")) {
            System.out.println(tridentTuple);
        }
    }
}




//output:
[5, 8187]
[1, 6381]
[2, 6090]
[4, 7731]
[3, 6131]
[7, 5988]
[6, 6400]

