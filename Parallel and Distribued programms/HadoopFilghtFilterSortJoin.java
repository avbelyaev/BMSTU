//hadoop lab3 (sort + join)
//задание:
//Требуется связать наборы данных по коду аэропорта прибытия: DEST_AEROPORT_ID
//Для каждого аэропорта требуется определить среднее, минимальное и максимальное время задержки для всех прибывающих рейсов.




//Запуск:
//hadoop fs -rmr output
//mvn package
//export HADOOP_CLASSPATH=target/hadoop-examples-1.0-SNAPSHOT.jar
//hadoop FlightSort 664600583_T_ONTIME_sample.csv L_AIRPORT_ID.csv output
//hadoop fs -copyToLocal output





//FlightSort.java
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.MultipleInputs;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.IntWritable;

public class FlightSort {
    public static void main(String[] args) throws Exception { 
        if (args.length != 3) {
            System.err.println("Usage: FlightSort-and-JoinApp <input path flight> <input path airport> <output path>");
            System.exit(-1); 
        } 

        Job job = Job.getInstance();

        job.setJarByClass(FlightSort.class);
        job.setJobName("Flight Sort-and-Join Job (hadlab3)");

        MultipleInputs.addInputPath(job, new Path(args[0]), TextInputFormat.class, FlightMapper.class);
        MultipleInputs.addInputPath(job, new Path(args[1]), TextInputFormat.class, AirportMapper.class);

        FileOutputFormat.setOutputPath(job, new Path(args[2]));

        job.setPartitionerClass(FlightPartitioner.class);
        job.setGroupingComparatorClass(FlightComparator.class);
        job.setReducerClass(FlightJoinReducer.class);
        job.setMapOutputKeyClass(FlightWritableComparable.class);
        job.setMapOutputValueClass(Text.class);

        //Mapper<.., .., KEYOUT, VALUEOUT>
        job.setOutputKeyClass(FlightWritableComparable.class);
        job.setOutputValueClass(Text.class);

        job.setNumReduceTasks(2);
        System.exit(job.waitForCompletion(true) ? 0 : 1); 
    } 
} 





//AirportMapper.java
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;

public class AirportMapper extends Mapper<LongWritable, Text, FlightWritableComparable, Text> {

    @Override
    protected void map(LongWritable key, Text value, Context context)
            throws IOException, InterruptedException {

        String[] columns = value.toString().replace("\"", "").split(",(?! )");

        if (!columns[0].equals("") &&
                !columns[0].equals("Code") &&
                !columns[1].equals("")) {



            FlightWritableComparable entry_writcom = new FlightWritableComparable();
            entry_writcom.setFlag(0);
            entry_writcom.setDest_airport_id(Integer.parseInt(columns[0]));

            Text airportname_text = new Text(columns[1]);

            context.write(entry_writcom, airportname_text);

        }

    }
}





//FlightMapper.java
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;

//Mapper<KEYIN, VALUEIN, KEYOUT, VALUEOUT>
public class FlightMapper extends Mapper<LongWritable, Text, FlightWritableComparable, Text> {

    @Override
    protected void map(LongWritable key, Text value, Context context)
            throws IOException, InterruptedException {

        String[] columns = value.toString().split(",");

        //flight
        if (!columns[0].equals("\"YEAR\"") && !columns[19].equals("")) {

            //not cancelled
            if ((float) 0 == Float.parseFloat(columns[19])) {

                //delayed
                if (!columns[18].equals("") &&
                        ((float)0 < Float.parseFloat(columns[18])) &&
                        !columns[14].equals("")) {



                    FlightWritableComparable entry_writcom = new FlightWritableComparable();
                    entry_writcom.setFlag(1);
                    entry_writcom.setDest_airport_id(Integer.parseInt(columns[14]));

                    Text delay_text = new Text(columns[18]);

                    context.write(entry_writcom, delay_text);

                }
            }
        }

    }
}






//FlightWritableComparable.java
import org.apache.hadoop.io.WritableComparable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class FlightWritableComparable implements WritableComparable {
    private int flag;
    private int dest_airport_id;

    public void setFlag(int flag) {
        this.flag = flag;
    }

    public void setDest_airport_id(int dest_airport_id) {
        this.dest_airport_id = dest_airport_id;
    }

    public int getDest_airport_id() {
        return dest_airport_id;
    }

    @Override
    public int compareTo(Object o) {

        FlightWritableComparable that = (FlightWritableComparable) o;

        int that_flag = that.flag;
        int that_dest_airport_id = that.dest_airport_id;

        if (this.dest_airport_id > that_dest_airport_id) {
            return 1;
        }
        if (this.dest_airport_id < that_dest_airport_id) {
            return -1;
        }

        if (this.flag > that_flag) {
            return 1;
        }
        if (this.flag < that_flag) {
            return -1;
        }

        return 0;

    }

    @Override
    public void write(DataOutput dataOutput) throws IOException {
        dataOutput.writeInt(flag);
        dataOutput.writeInt(dest_airport_id);
    }

    @Override
    public void readFields(DataInput dataInput) throws IOException {
        flag = dataInput.readInt();
        dest_airport_id = dataInput.readInt();
    }

    @Override
    public int hashCode() {

        int hash = this.toString().hashCode();

        if (hash < 0) {
            return -hash;
        }

        return  hash;
    }

    @Override
    public String toString() {
        return "WritComp{" +
                "flag=" + flag +
                ", dest_airport_id=" + dest_airport_id +
                '}';
    }

    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FlightWritableComparable that = (FlightWritableComparable) o;

        return flag == that.flag && dest_airport_id == that.dest_airport_id;

    }
}






//FlightPartitioner.java
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Partitioner;

public class FlightPartitioner extends Partitioner<FlightWritableComparable, Text> {

    @Override
    public int getPartition(FlightWritableComparable key, Text value, int numPartitions) {
        return key.getDest_airport_id() % numPartitions;
    }
}





//FlightComparator.java
import org.apache.hadoop.io.WritableComparable;
import org.apache.hadoop.io.WritableComparator;

public class FlightComparator extends WritableComparator {

    public FlightComparator() {
        super(FlightWritableComparable.class, true);
    }

    @Override
    public int compare(WritableComparable a1, WritableComparable b1) {

        FlightWritableComparable a = (FlightWritableComparable) a1;
        FlightWritableComparable b = (FlightWritableComparable) b1;

        return Integer.compare(a.getDest_airport_id(), b.getDest_airport_id());
    }
}






//FlightReducer.java
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import java.io.IOException;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//Reducer<KEYIN, VALUEIN, KEYOUT, VALUEOUT>
public class FlightJoinReducer extends Reducer<FlightWritableComparable, Text, String, String> {

    @Override
    protected void reduce(FlightWritableComparable key, Iterable<Text> values, Context context)
            throws IOException, InterruptedException {

        String airport_name = "AirPort:";
        String airport_stats = "\tDelays:";
        //String tmp = "";

        Iterator<Text> val_iterator = values.iterator();

        airport_name += val_iterator.next().toString();

        if (val_iterator.hasNext()) {

            int i = 0;
            float tmp = Float.parseFloat(val_iterator.next().toString());
            float sum = (float) 0, min = tmp, max = tmp;

            while (val_iterator.hasNext()) {

                float curr_val = Float.parseFloat(val_iterator.next().toString());

                if (curr_val > max) {
                    max = curr_val;
                }
                if (curr_val < min) {
                    min = curr_val;
                }
                sum += curr_val;
                i++;

            }

            airport_stats += "Min:" + min;
            airport_stats += ", Max:" + max;
            airport_stats += ", Avg:" + (sum / (float)i);

            context.write(airport_name, airport_stats);

        }
    }
}






//Example:
AirPort:Abilene, TX: Abilene Regional		Delays:Min:1.0, Max:76.0, Avg:19.666666
AirPort:Albuquerque, NM: Albuquerque International Sunport		Delays:Min:1.0, Max:126.0, Avg:19.852942
AirPort:Albany, GA: Southwest Georgia Regional		Delays:Min:12.0, Max:34.0, Avg:21.0
AirPort:Atlantic City, NJ: Atlantic City International		Delays:Min:2.0, Max:147.0, Avg:58.153847
AirPort:Kodiak, AK: Kodiak Airport		Delays:Min:54.0, Max:76.0, Avg:76.0
AirPort:Augusta, GA: Augusta Regional at Bush Field		Delays:Min:1.0, Max:246.0, Avg:50.8
AirPort:Waterloo, IA: Waterloo Regional		Delays:Min:2.0, Max:223.0, Avg:60.0
AirPort:Aspen, CO: Aspen Pitkin County Sardy Field		Delays:Min:5.0, Max:136.0, Avg:36.264706
AirPort:Appleton, WI: Outagamie County Regional		Delays:Min:3.0, Max:146.0, Avg:54.333332
AirPort:Scranton/Wilkes-Barre, PA: Wilkes Barre Scranton International		Delays:Min:1.0, Max:20.0, Avg:10.666667
AirPort:Billings, MT: Billings Logan International		Delays:Min:5.0, Max:58.0, Avg:31.142857
AirPort:Bellingham, WA: Bellingham International		Delays:Min:3.0, Max:20.0, Avg:3.0
AirPort:Beaumont/Port Arthur, TX: Jack Brooks Regional		Delays:Min:14.0, Max:36.0, Avg:21.0
AirPort:Aguadilla, PR: Rafael Hernandez		Delays:Min:3.0, Max:109.0, Avg:7.0
AirPort:Barrow, AK: Wiley Post/Will Rogers Memorial		Delays:Min:2.0, Max:306.0, Avg:67.333336
AirPort:Buffalo, NY: Buffalo Niagara International		Delays:Min:1.0, Max:128.0, Avg:23.833334
AirPort:Burbank, CA: Bob Hope		Delays:Min:1.0, Max:192.0, Avg:25.863636
AirPort:Columbia, SC: Columbia Metropolitan		Delays:Min:1.0, Max:119.0, Avg:30.185184
AirPort:Akron, OH: Akron-Canton Regional		Delays:Min:1.0, Max:83.0, Avg:20.6
AirPort:Cedar City, UT: Cedar City Regional		Delays:Min:1.0, Max:15.0, Avg:1.0
AirPort:Cordova, AK: Merle K Mudhole Smith		Delays:Min:1.0, Max:151.0, Avg:1.0

