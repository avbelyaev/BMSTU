//hadoop lab2
//задание:
//а.Требуется отфильтровать рейсы которые были отменены или прошли с опозданием
//б. Данные, прошедшие фильтр, требуется отсортировать по времени опоздания, аэропорту назначения и времени перелета


//запуск:
//hadoop fs -rmr output
//mvn package
//export HADOOP_CLASSPATH=target/hadoop-examples-1.0-SNAPSHOT.jar
//hadoop FlightSort 664600583_T_ONTIME_sample.csv output
//hadoop fs -copyToLocal output



//FlightSort.java
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.IntWritable;

public class FlightSort {
    public static void main(String[] args) throws Exception { 
        if (args.length != 2) { 
            System.err.println("Usage: FlightSortApp <input path> <output path>");
            System.exit(-1); 
        } 

        Job job = Job.getInstance();

        job.setJarByClass(FlightSort.class);
        job.setJobName("Flight Sort Job (hadlab2)");

        FileInputFormat.addInputPath(job, new Path(args[0]));
        FileOutputFormat.setOutputPath(job, new Path(args[1]));

        job.setMapperClass(FlightMapper.class);
        job.setPartitionerClass(FlightPartitioner.class);
        job.setReducerClass(FlightReducer.class);
        job.setGroupingComparatorClass(FlightComparator.class);

        //Mapper<.., .., KEYOUT, VALUEOUT>
        job.setOutputKeyClass(FlightWritableComparable.class);
        job.setOutputValueClass(Text.class);

        job.setNumReduceTasks(2);
        System.exit(job.waitForCompletion(true) ? 0 : 1); 
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

        //flight (not metadata)
        if (!columns[0].equals("\"YEAR\"") && !columns[19].equals("")) {

            //not cancelled
            if ((float)0 == Float.parseFloat(columns[19])) {

                //delayed
                if (!columns[18].equals("") && ((float)0 != Float.parseFloat(columns[18])) &&
                        !columns[14].equals("") &&
                        !columns[21].equals("")) {

                    FlightWritableComparable entry = new FlightWritableComparable();

                    entry.setArr_delay(Float.parseFloat(columns[18]));
                    entry.setDest_airport_id(Integer.parseInt(columns[14]));
                    entry.setAir_time(Float.parseFloat(columns[21]));
                    entry.setCancelled((float)0);

                    context.write(entry, value);

                }

            //cancelled
            } else {

                //so it has no ArrDelay, no AirTime, no DestAirport, only OriginAriport
                if (!columns[11].equals("")) {

                    FlightWritableComparable entry = new FlightWritableComparable();

                    entry.setArr_delay((float)0);
                    entry.setDest_airport_id(Integer.parseInt(columns[11]));
                    entry.setAir_time((float)0);
                    entry.setCancelled((float)1);

                    context.write(entry, value);

                }

            }

        }

    }
}




//FlightReducer.java
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import java.io.IOException;
import java.util.Iterator;

//Reducer<KEYIN, VALUEIN, KEYOUT, VALUEOUT>
public class FlightReducer extends Reducer<FlightWritableComparable, Text, String, Text> {
    @Override
    protected void reduce(FlightWritableComparable key, Iterable<Text> values, Context context)
            throws IOException, InterruptedException {

        for (Text t : values){

            String brief_info = "Flight: Brief: [" + key.getBriefInfo() + "]";

            String tmp = t.toString();
            Text full_info = new Text("Full: [" + tmp + "]");

            context.write(brief_info, full_info);
        }

    }
}




//FlightWritableComparable.java
import org.apache.hadoop.io.WritableComparable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class FlightWritableComparable implements WritableComparable {

    private float arr_delay;
    private int dest_airport_id;
    private float air_time;
    private float cancelled;

    public void setArr_delay(float arr_delay) {
        this.arr_delay = arr_delay;
    }

    public void setDest_airport_id(int dest_airport_id) {
        this.dest_airport_id = dest_airport_id;
    }

    public void setAir_time(float air_time) {
        this.air_time = air_time;
    }

    public void setCancelled(float cancelled) {
        this.cancelled = cancelled;
    }

    //for comparator's compare
    public int compareDelayTime(Object o) {
        FlightWritableComparable that = (FlightWritableComparable) o;

        float that_arr_delay = that.arr_delay;

        if (this.arr_delay > that_arr_delay) {
            return 1;
        }
        if (this.arr_delay < that_arr_delay) {
            return -1;
        }

        return 0;
    }

    @Override
    public int compareTo(Object o) {

        FlightWritableComparable that = (FlightWritableComparable) o;

        float that_arr_delay = that.arr_delay;
        int that_dest_airport_id = that.dest_airport_id;
        float that_air_time = that.air_time;
        float that_cancelled = that.cancelled;


        if (this.cancelled > that_cancelled) {
            return 1;
        }
        if (this.cancelled < that_cancelled) {
            return -1;
        }

        if (this.arr_delay > that_arr_delay) {
            return 1;
        }
        if (this.arr_delay < that_arr_delay) {
            return -1;
        }

        if (this.dest_airport_id > that_dest_airport_id) {
            return 1;
        }
        if (this.dest_airport_id < that_dest_airport_id) {
            return -1;
        }

        if (this.air_time > that_air_time) {
            return 1;
        }
        if (this.air_time < that_air_time) {
            return -1;
        }

        return 0;

    }

    @Override
    public void write(DataOutput dataOutput) throws IOException {
        dataOutput.writeFloat(arr_delay);
        dataOutput.writeInt(dest_airport_id);
        dataOutput.writeFloat(air_time);
        dataOutput.writeFloat(cancelled);
    }

    @Override
    public void readFields(DataInput dataInput) throws IOException {
        arr_delay = dataInput.readFloat();
        dest_airport_id = dataInput.readInt();
        air_time = dataInput.readFloat();
        cancelled = dataInput.readFloat();
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
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FlightWritableComparable that = (FlightWritableComparable) o;

        if (Float.compare(that.arr_delay, arr_delay) != 0) return false;
        if (dest_airport_id != that.dest_airport_id) return false;
        if (Float.compare(that.air_time, air_time) != 0) return false;
        return Float.compare(that.cancelled, cancelled) == 0;

    }

    //just in case
    public String getBriefInfo() {
        return "Delay:" + arr_delay +
                ", AirportID:" + dest_airport_id +
                ", AirTime:" + air_time +
                ", Cancelled:" + cancelled;
    }

    @Override
    public String toString() {
        return "FlightWritableComparable{" +
                "arr_delay=" + arr_delay +
                ", dest_airport_id=" + dest_airport_id +
                ", air_time=" + air_time +
                ", cancelled=" + cancelled +
                '}';
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

        return a.compareDelayTime(b);
    }
}




//FlightPartitioner.java
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Partitioner;

public class FlightPartitioner extends Partitioner<FlightWritableComparable, Text> {

    @Override
    public int getPartition(FlightWritableComparable key, Text value, int numPartitions) {
        return key.hashCode() % numPartitions;
    }
}





//output example:
Flight: Brief: [Delay:1.0, AirportID:10140, AirTime:78.0, Cancelled:0.0]	Full: [2015,1,1,10,6,2015-01-10,"OO",20304,"OO","N170PQ","4499",14869,1486903,34614,10140,"1844","1853",1.00,1.00,0.00,"",78.00,493.00,]
Flight: Brief: [Delay:1.0, AirportID:10140, AirTime:85.0, Cancelled:0.0]	Full: [2015,1,1,25,7,2015-01-25,"WN",19393,"WN","N7739A","4285",14679,1467903,33570,10140,"2317","2321",1.00,1.00,0.00,"",85.00,628.00,]
Flight: Brief: [Delay:1.0, AirportID:10140, AirTime:153.0, Cancelled:0.0]	Full: [2015,1,1,10,6,2015-01-10,"WN",19393,"WN","N224WN","2739",14747,1474703,30559,10140,"1530","1541",1.00,1.00,0.00,"",153.00,1180.00,]
Flight: Brief: [Delay:1.0, AirportID:10279, AirTime:80.0, Cancelled:0.0]	Full: [2015,1,1,5,1,2015-01-05,"EV",20366,"EV","N14514","5792",12266,1226603,31453,10279,"1607","1613",1.00,1.00,0.00,"",80.00,517.00,]
Flight: Brief: [Delay:1.0, AirportID:10299, AirTime:198.0, Cancelled:0.0]	Full: [2015,1,1,13,2,2015-01-13,"AS",19930,"AS","N303AS","85",14747,1474703,30559,10299,"1307","1312",1.00,1.00,0.00,"",198.00,1448.00,]
Flight: Brief: [Delay:1.0, AirportID:10397, AirTime:45.0, Cancelled:0.0]	Full: [2015,1,1,23,5,2015-01-23,"DL",19790,"DL","N978DL","1956",11481,1148102,31481,10397,"1915","1927",1.00,1.00,0.00,"",45.00,240.00,]
Flight: Brief: [Delay:1.0, AirportID:10397, AirTime:47.0, Cancelled:0.0]	Full: [2015,1,1,12,1,2015-01-12,"DL",19790,"DL","N986DL","2422",11057,1105703,31057,10397,"1044","1054",1.00,1.00,0.00,"",47.00,226.00,]
Flight: Brief: [Delay:1.0, AirportID:10397, AirTime:56.0, Cancelled:0.0]	Full: [2015,1,1,30,5,2015-01-30,"EV",20366,"EV","N861AS","4962",14574,1457402,34574,10397,"1524","1530",1.00,1.00,0.00,"",56.00,357.00,]
Flight: Brief: [Delay:1.0, AirportID:10397, AirTime:63.0, Cancelled:0.0]	Full: [2015,1,1,4,7,2015-01-04,"WN",19393,"WN","N947WN","888",13495,1349503,33495,10397,"0853","0901",1.00,1.00,0.00,"",63.00,425.00,]
Flight: Brief: [Delay:1.0, AirportID:10397, AirTime:69.0, Cancelled:0.0]	Full: [2015,1,1,9,5,2015-01-09,"DL",19790,"DL","N357NW","755",11066,1106603,31066,10397,"1636","1641",1.00,1.00,0.00,"",69.00,447.00,]
Flight: Brief: [Delay:1.0, AirportID:10397, AirTime:70.0, Cancelled:0.0]	Full: [2015,1,1,16,5,2015-01-16,"DL",19790,"DL","N907DE","2313",14122,1412202,30198,10397,"2116","2122",1.00,1.00,0.00,"",70.00,526.00,]

