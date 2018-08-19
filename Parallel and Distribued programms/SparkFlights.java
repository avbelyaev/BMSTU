//hadlab5: Spark
//задание:
//Требуется определить для пары <аэропорт отлета, аэропорт прибытия> максимальное время опоздания, процент опоздавших+отмененных рейсов.
//Также требуется связать полученную таблицу с названиями аэропортов.




//запуск:
//spark-submit --class SparkLab --master yarn-client --num-executors 3 target/spark-examples-1.0-SNAPSHOT.jar




//SparkLab.java
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import scala.Tuple2;

import java.io.Serializable;
import java.lang.*;
import java.lang.Float;
import java.lang.Integer;
import java.lang.Object;
import java.util.Arrays;
import java.lang.Boolean;
import java.lang.Exception;
import java.lang.Iterable;
import java.lang.Override;
import java.lang.String;
import java.lang.System;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Pattern;

public class SparkLab {

    public static class FlightSerializable implements Serializable {
        //just in case
        public int origin_airport_id;
        public int dest_airport_id;
        //real data
        public float cancelled;
        public float arr_delay_new;

        public FlightSerializable(int origin_airport_id, int dest_airport_id, float cancelled, float arr_delay_new) {
            this.origin_airport_id = origin_airport_id;
            this.dest_airport_id = dest_airport_id;
            this.cancelled = cancelled;
            this.arr_delay_new = arr_delay_new;
        }

        public FlightSerializable(int origin_airport_id, int dest_airport_id, float cancelled) {
            this.origin_airport_id = origin_airport_id;
            this.dest_airport_id = dest_airport_id;
            this.cancelled = cancelled;
        }
    }

    public static void main(String[] args) {

        SparkConf conf = new SparkConf().setAppName("sparklab5");
        JavaSparkContext sc = new JavaSparkContext(conf);

        JavaRDD<String> flights = sc.textFile("664600583_T_ONTIME_sample.csv");
        JavaRDD<String> airports = sc.textFile("L_AIRPORT_ID.csv");



	//load airport data from csv to rdd of strings
        JavaRDD<String> a_no_meta = airports.map(new Function<String, String>() {
            @Override
            public String call(String line) throws Exception {
                return line.replace("\"", "");
            }
        }).filter(new Function<String, java.lang.Boolean>() {
            @Override
            public Boolean call(String line2) throws Exception {
                return !line2.contains("Code") && !line2.contains("Description");
            }
        });



	//generate pairs: <airport_code, description> 
        JavaPairRDD<Integer, String> a_paired = a_no_meta.mapToPair(new PairFunction<String, Integer, String>() {
            @Override
            public Tuple2<Integer, String> call(String line) throws Exception {

                String[] columns = line.split(",(?! )");
		//dummy data
                int airport_code = -1;
                String description = "dummy";

                if (!columns[0].equals("") && !columns[1].equals("")) {
                    airport_code = Integer.parseInt(columns[0]);
                    description = columns[1];
                }
                return new Tuple2<Integer, String>(airport_code, description);
            }
        });



	//get same data as Map and broadcast
        final Map<Integer, String> a_map = a_paired.collectAsMap();

        final Broadcast<Map<Integer, String>> airports_broadcasted = sc.broadcast(a_map);






	//load flight data from csv to rdd of strings
        JavaRDD<String> f_no_meta = flights.map(new Function<String, String>() {
            @Override
            public String call(String line) throws Exception {
                return line.replace("\"", "");
            }
        }).filter(new Function<String, Boolean>() {
            @Override
            public Boolean call(String line2) throws Exception {
                return !line2.contains("YEAR");
            }
        });


	//generate corteges of 
	//		KEY: cortege of: <origin_airport_id, dest_airport_id>
	//		VALUE: <FlightSerializable> container for data (airports, delay, cancellation_mark) as Serializable object
        JavaPairRDD<Tuple2<Integer, Integer>, FlightSerializable> f_paired = f_no_meta.mapToPair(new PairFunction<String, Tuple2<Integer, Integer>, FlightSerializable>() {
            @Override
            public Tuple2<Tuple2<Integer, Integer>, FlightSerializable> call(String line) throws Exception {
                String[] columns = line.split(",");

                //dummy data to return
                int origin_ap = -1;
                int dest_ap = -1;
                float cnc = (float)1;

                FlightSerializable flight_data_storage =
                        new FlightSerializable(origin_ap, dest_ap, cnc);


                //have origin_ap, dest_ap, cancellation_flag
                if (!columns[11].equals("") &&
                        !columns[14].equals("") &&
                        !columns[19].equals("")) {

                    origin_ap = Integer.parseInt(columns[11]);
                    dest_ap = Integer.parseInt(columns[14]);

                    flight_data_storage.origin_airport_id = Integer.parseInt(columns[11]);
                    flight_data_storage.dest_airport_id = Integer.parseInt(columns[14]);


                    //not cancelled
                    if ((float)0 == Float.parseFloat(columns[19]) &&
                            !columns[18].equals("")) {

                        flight_data_storage.cancelled = (float)0;
                        flight_data_storage.arr_delay_new = Float.parseFloat(columns[18]);

                    }
                    //cancelled
                    if ((float)1 == Float.parseFloat(columns[19])) {

                        flight_data_storage.cancelled = (float)1;
                        //777 = dummy
                        flight_data_storage.arr_delay_new = (float)777;

                    }
                }

                return new Tuple2<Tuple2<Integer, Integer>, FlightSerializable>
                        (new Tuple2<Integer, Integer>(origin_ap, dest_ap),
                                flight_data_storage);

            }
        });

	//pairs example:
	//((12478,12892),SparkLab$FlightSerializable@7aec3716)
	//((12478,12892),SparkLab$FlightSerializable@47d88d7b)
	//((12478,12892),SparkLab$FlightSerializable@5a30ab16)




        JavaPairRDD<Tuple2<Integer, Integer>, Iterable<FlightSerializable>> f_groupped =
                f_paired.groupByKey();





	//simplify iterable of Serialized objects to string in same cortage:
	//		KEY: cortege of: <origin_airport_id, dest_airport_id>
	//		VALUE: String based on data from all serialized objects that were groupped by that KEY
        JavaPairRDD<Tuple2<Integer, Integer>, String> f_mapped = f_groupped.mapToPair(
                new PairFunction<Tuple2<Tuple2<Integer, Integer>, Iterable<FlightSerializable>>, Tuple2<Integer, Integer>, String>() {
            @Override
            public Tuple2<Tuple2<Integer, Integer>, String> call(Tuple2<Tuple2<Integer, Integer>, Iterable<FlightSerializable>> value) throws Exception {

                Tuple2<Integer, Integer> airport_bundle = value._1;
                Iterable<FlightSerializable> data_containers_bundle = value._2;

                String ret = "";
                int flights_quantity = 0;


                int num_of_positive_delays = 0;
                int num_of_cancellations = 0;
                float max_delay = (float)0;


                for (FlightSerializable t : data_containers_bundle) {

                    if ((float)0 == t.cancelled) {

                        if ((float)0 < t.arr_delay_new) {
                            if (t.arr_delay_new > max_delay) {
                                max_delay = t.arr_delay_new;
                            }
                            num_of_positive_delays++;
                        }
                    }

                    if ((float)1 == t.cancelled) {
                        num_of_cancellations++;
                    }

                    flights_quantity++;
                }

                ret += "Flights Qty:" + flights_quantity + ", Stats:";

                float delay_percentage = (float)num_of_positive_delays / (float)flights_quantity * (float)100;
                float cancellation_percentage = (float)num_of_cancellations / (float)flights_quantity * (float)100;

                ret += "Cancelled:" + num_of_cancellations + "/" + flights_quantity + "[" + cancellation_percentage + "%], ";

                if ((float)0 == max_delay) {
                    ret += "Delayed: - [0.0%], Max_Delay: - [0.0]";
                } else {
                    ret += "Delayed:" + num_of_positive_delays + "/" + flights_quantity + "[" + delay_percentage + "%], ";
                    ret += "Max_Delay:" + max_delay;
                }

                return new Tuple2<Tuple2<Integer, Integer>, String>(airport_bundle, ret);
            }
        });

	//example:
	//((10299,10926),Flights Qty:8, Stats:Cancelled:0/8[0.0%], Delayed:2/8[25.0%], Max_Delay:151.0)
	//((14107,13891),Flights Qty:21, Stats:Cancelled:1/21[4.7619047%], Delayed:9/21[42.857143%], Max_Delay:67.0)
	//((12266,10721),Flights Qty:1, Stats:Cancelled:0/1[0.0%], Delayed: - [0.0%], Max_Delay: - [0.0])



	//get airport names from previous broadcasted Map and concatenate all data to same string
        JavaRDD<String> final_flight_dataset = f_mapped.map(new Function<Tuple2<Tuple2<Integer, Integer>, String>, String>() {
            @Override
            public String call(Tuple2<Tuple2<Integer, Integer>, String> value) throws Exception {

                Map<Integer, String> a_names_broadcasted_map = airports_broadcasted.value();
                Tuple2<Integer, Integer> airport_code_bundle = value._1;
                String airport_bundle_stats = value._2;

                String origin_ap_name = a_names_broadcasted_map.get(airport_code_bundle._1);
                String dest_ap_name = a_names_broadcasted_map.get(airport_code_bundle._2);

                return "F: " + origin_ap_name + " ->> " + dest_ap_name + ", \t" + airport_bundle_stats;

            }
        });

        final_flight_dataset.saveAsTextFile("hdfs://localhost:9000/user/anthony/spar");



        System.out.println("\nSUCCESSFULLY EXECUTED:" +
                ", in(flights):" + f_paired.count() +
                ", out(bundles):" + final_flight_dataset.count() + "\n");
    }

}





//Output example:
F: Honolulu, HI: Honolulu International ->> Dallas/Fort Worth, TX: Dallas/Fort Worth International, 	Flights Qty:7, Stats:Cancelled:0/7[0.0%], Delayed:5/7[71.42857%], Max_Delay:70.0
F: Midland/Odessa, TX: Midland International ->> Dallas, TX: Dallas Love Field, 	Flights Qty:11, Stats:Cancelled:0/11[0.0%], Delayed:4/11[36.363636%], Max_Delay:29.0
F: Orlando, FL: Orlando International ->> Chicago, IL: Chicago O'Hare International, 	Flights Qty:41, Stats:Cancelled:0/41[0.0%], Delayed:18/41[43.90244%], Max_Delay:94.0
F: Detroit, MI: Detroit Metro Wayne County ->> Chicago, IL: Chicago Midway International, 	Flights Qty:23, Stats:Cancelled:1/23[4.347826%], Delayed:12/23[52.173912%], Max_Delay:78.0
F: Denver, CO: Denver International ->> San Antonio, TX: San Antonio International, 	Flights Qty:10, Stats:Cancelled:0/10[0.0%], Delayed:7/10[70.0%], Max_Delay:167.0
F: Charlotte, NC: Charlotte Douglas International ->> Phoenix, AZ: Phoenix Sky Harbor International, 	Flights Qty:30, Stats:Cancelled:0/30[0.0%], Delayed:14/30[46.666668%], Max_Delay:220.0
F: Chicago, IL: Chicago O'Hare International ->> Baltimore, MD: Baltimore/Washington International Thurgood Marshall, 	Flights Qty:22, Stats:Cancelled:1/22[4.5454545%], Delayed:13/22[59.090908%], Max_Delay:475.0
F: Las Vegas, NV: McCarran International ->> Wichita, KS: Wichita Dwight D Eisenhower National, 	Flights Qty:4, Stats:Cancelled:0/4[0.0%], Delayed:1/4[25.0%], Max_Delay:31.0
F: Chicago, IL: Chicago Midway International ->> Los Angeles, CA: Los Angeles International, 	Flights Qty:13, Stats:Cancelled:0/13[0.0%], Delayed:6/13[46.153847%], Max_Delay:93.0
F: Miami, FL: Miami International ->> Washington, DC: Ronald Reagan Washington National, 	Flights Qty:28, Stats:Cancelled:1/28[3.5714288%], Delayed:12/28[42.857143%], Max_Delay:44.0
F: San Francisco, CA: San Francisco International ->> Kona, HI: Kona International Airport at Keahole, 	Flights Qty:8, Stats:Cancelled:0/8[0.0%], Delayed:3/8[37.5%], Max_Delay:28.0

