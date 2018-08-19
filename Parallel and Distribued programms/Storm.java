//hadlab7: Storm
Лабораторная работа 7. Разработка простой топологии storm.
Задача:
-Требуется разработать топологию для составления частотного словаря файлов.
-Требуется разработать Spout который будет опрашивать директорию и в случае
	обнаружения в ней новых файлов читать ее построчно и генерировать tuple в
	исходящий поток words для каждой строки. После завершения чтения файла
	Spout должен выдать tuple в исходящий поток sync. Также после окончания
	файла требуется перенести файл в папку для обработанных файлов
-Требуется разработать Bolt Splitter который будет принимать строку и разбивать ее на слова.
-Требуется разработать Bolt Counter который будет принимать слова из
	входящего потока Splitter и вести частотный словарь. Также по команде от
	входного потока sync — печатать текущий словарь на экране и обнулять его.



//запуск:
//mvn package
//mvn exec:java -Dexec.mainClass="StormLab"



//StormLab.java
import backtype.storm.Config;
import backtype.storm.LocalCluster;
import backtype.storm.topology.TopologyBuilder;
import backtype.storm.tuple.Fields;

public class StormLab {

    public static void main(String[] args) throws Exception {

        TopologyBuilder builder = new TopologyBuilder();
        builder.setSpout("generator", new PollSpout());

        builder.setBolt("splitter", new BoltSplitter(), 10)
                .shuffleGrouping("generator", "words");

        builder.setBolt("counter", new BoltCounter(), 1)
                .fieldsGrouping("splitter", new Fields("word"))
                .allGrouping("generator", "sync");

        Config config = new Config();
        config.setDebug(false);

        LocalCluster cluster = new LocalCluster();
        cluster.submitTopology("Frequency Dictionary", config, builder.createTopology());
    }
}






//PollSpout.java
import backtype.storm.spout.SpoutOutputCollector;
import backtype.storm.task.TopologyContext;
import backtype.storm.topology.OutputFieldsDeclarer;
import backtype.storm.topology.base.BaseRichSpout;
import backtype.storm.tuple.Fields;
import backtype.storm.tuple.Values;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import org.apache.storm.shade.com.google.common.base.Charsets;

import java.io.*;
import java.util.Map;

public class PollSpout extends BaseRichSpout {
    private SpoutOutputCollector output;
    private boolean reading_currently;
    private File dir, current_file;
    private BufferedReader reader;
    private int ack_num, emit_num;

    public PollSpout (){
        emit_num = 0;
        ack_num = 0;
        reading_currently = false;
        dir = new File("/home/anthony/hadlab7/in_files");

		//should exist 2 directories: 
			//in_files with input file in it
			//out_files will contain input file after its processed

    }

    public void open(Map map, TopologyContext topologyContext, SpoutOutputCollector spoutOutputCollector) {
        output = spoutOutputCollector;
    }

    public void declareOutputFields(OutputFieldsDeclarer outputFieldsDeclarer) {
        outputFieldsDeclarer.declareStream("words", new Fields("words"));
        outputFieldsDeclarer.declareStream("sync", new Fields());
    }

    public void nextTuple() {

        if (reading_currently) {
            try {

                String line = reader.readLine();
                if (null != line) {

                    output.emit("words", new Values(line), emit_num);
                    emit_num++;

                } else {
                    if(ack_num == emit_num) {

                        File dest_file = new File("/home/anthony/hadlab7/out_files/" + current_file.getName());
                        Files.move(current_file, dest_file);

                        output.emit("sync", Lists.newArrayList());
                        ack_num = 0;
                        emit_num = 0;
                        reading_currently = false;

                    }
                }
            } catch (IOException e) { e.printStackTrace(); }
        } else {

            File files_list[] = dir.listFiles();
            if (null == files_list || 0 == files_list.length) {

                backtype.storm.utils.Utils.sleep(100);

            } else {
                try {

                    current_file = files_list[0];
                    reader = new BufferedReader(new InputStreamReader(
                            new FileInputStream(current_file), Charsets.UTF_8));

                } catch (FileNotFoundException e) { e.printStackTrace(); }
                reading_currently = true;
            }
        }
    }

    public void ack(Object msgId) {
        ack_num++;
    }
}





//BoltSplitter.java
import backtype.storm.task.OutputCollector;
import backtype.storm.task.TopologyContext;
import backtype.storm.topology.OutputFieldsDeclarer;
import backtype.storm.topology.base.BaseRichBolt;
import backtype.storm.tuple.Fields;
import backtype.storm.tuple.Tuple;
import com.google.common.collect.Lists;

import java.util.Map;

public class BoltSplitter extends BaseRichBolt {
    private OutputCollector collector;

    public void prepare(Map stormConf, TopologyContext context, OutputCollector collector) {
        this.collector = collector;
    }

    public void declareOutputFields(OutputFieldsDeclarer declarer) {
        declarer.declare(new Fields("word"));
    }

    public void execute(Tuple input) {
        String words[] = input.getStringByField("words").split("[^a-zA-Zа-яА-Я]+");
        for (String w : words){
            collector.emit(input, Lists.newArrayList((Object) w));
        }
        collector.ack(input);
    }
}





//BoltCounter.java
import backtype.storm.task.OutputCollector;
import backtype.storm.task.TopologyContext;
import backtype.storm.topology.OutputFieldsDeclarer;
import backtype.storm.topology.base.BaseRichBolt;
import backtype.storm.tuple.Tuple;
import java.util.HashMap;
import java.util.Map;

public class BoltCounter extends BaseRichBolt {
    private OutputCollector output;
    private Map<String, Integer> table;

    public BoltCounter(){}

    public void declareOutputFields(OutputFieldsDeclarer declarer) {}

    public void prepare(Map stormConf, TopologyContext context, OutputCollector collector) {
        output = collector;
        table = new HashMap<String, Integer>();
    }

    public void execute(Tuple tuple) {
        try {
            if (tuple.getSourceStreamId().equals("sync")) {

                for (Map.Entry<String, Integer> s : table.entrySet()) {
                    System.out.println(s.getKey() + " : " + s.getValue());
                }
                table = new HashMap<String, Integer>();

            } else {

                String word = tuple.getStringByField("word");
                Integer counts = table.get(word);

                if (null == counts)
                    counts = 0;

                counts++;
                table.put(word, counts);
                output.ack(tuple);

            }
        } catch (Exception e) { e.printStackTrace(); }
    }
}

