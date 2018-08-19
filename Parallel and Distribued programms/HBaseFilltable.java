import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.TableName;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class FillTable {
    //col family
    private static final byte[] CF = "data".getBytes();
    //cols
    private static final byte[] ATTR0_YEAR = "year".getBytes();
    private static final byte[] ATTR1_QUARTER = "quarter".getBytes();
    private static final byte[] ATTR2_MONTH = "month".getBytes();
    private static final byte[] ATTR3_DAY_OF_MONTH = "day_of_month".getBytes();
    private static final byte[] ATTR4_DAY_OF_WEEK = "day_of_week".getBytes();
    private static final byte[] ATTR5_FL_DATE = "fl_date".getBytes();
    private static final int FL_DATE = 5;
    private static final byte[] ATTR6_UNIQUE_CARRIER = "unique_carrier".getBytes();
    private static final byte[] ATTR7_AIRLINE_ID = "airline_id".getBytes();
    private static final int AIRLINE_ID = 7;
    private static final byte[] ATTR8_CARRIER = "carrier".getBytes();
    private static final byte[] ATTR9_TAIL_NUM = "tail_num".getBytes();
    private static final byte[] ATTR10_FL_NUM = "fl_num".getBytes();
    private static final byte[] ATTR11_ORIGIN_AIRPORT_ID = "origin_airport_id".getBytes();
    private static final byte[] ATTR12_ORIGIN_AIRPORT_SEQ_ID = "origin_airport_seq_id".getBytes();
    private static final byte[] ATTR13_ORIGIN_CITY_MARKET_ID = "origin_city_market_id".getBytes();
    private static final byte[] ATTR14_DEST_AIRPORT_ID = "dest_airport_id".getBytes();
    private static final byte[] ATTR15_WHEELS_ON = "wheels_on".getBytes();
    private static final byte[] ATTR16_ARR_TIME = "arr_time".getBytes();
    private static final byte[] ATTR17_ARR_DELAY = "arr_delay".getBytes();
    private static final byte[] ATTR18_ARR_DELAY_NEW = "arr_delay_new".getBytes();
    private static final byte[] ATTR19_CANCELLED = "cancelled".getBytes();
    private static final byte[] ATTR20_CANCELLATION_CODE = "cancellation_code".getBytes();
    private static final byte[] ATTR21_AIR_TIME = "air_time".getBytes();
    private static final byte[] ATTR22_DISTANCE = "distance".getBytes();


    public static void main(String[] args) throws IOException {

        Configuration config = HBaseConfiguration.create();
        config.set("hbase.zookeper.quorum", "localhost");

        //hbase> create "flights", "data"
        Connection connection = ConnectionFactory.createConnection(config);
        Table table = connection.getTable(TableName.valueOf("flights"));

        String data_path = "/home/anthony/hadlab6/664600583_T_ONTIME_sample.csv";
        BufferedReader reader = new BufferedReader(new FileReader(data_path));

        int row_num = 0;
        while (true) {

            String line = reader.readLine();
            if (null == line) {
                break;
            }

            String[] columns = line.replace("\"", "").split(",");

            if (!columns[0].equals("YEAR")) {

                Put put = new Put(Bytes.toBytes(
                        columns[FL_DATE] + "_" + columns[AIRLINE_ID] + "_" + row_num));

                put.addColumn(CF, ATTR0_YEAR, Bytes.toBytes(columns[0]));
                put.addColumn(CF, ATTR1_QUARTER, Bytes.toBytes(columns[1]));
                put.addColumn(CF, ATTR2_MONTH, Bytes.toBytes(columns[2]));
                put.addColumn(CF, ATTR3_DAY_OF_MONTH, Bytes.toBytes(columns[3]));
                put.addColumn(CF, ATTR4_DAY_OF_WEEK, Bytes.toBytes(columns[4]));
                put.addColumn(CF, ATTR5_FL_DATE, Bytes.toBytes(columns[5]));
                put.addColumn(CF, ATTR6_UNIQUE_CARRIER, Bytes.toBytes(columns[6]));
                put.addColumn(CF, ATTR7_AIRLINE_ID, Bytes.toBytes(columns[7]));
                put.addColumn(CF, ATTR8_CARRIER, Bytes.toBytes(columns[8]));
                put.addColumn(CF, ATTR9_TAIL_NUM, Bytes.toBytes(columns[9]));
                put.addColumn(CF, ATTR10_FL_NUM, Bytes.toBytes(columns[10]));
                put.addColumn(CF, ATTR11_ORIGIN_AIRPORT_ID, Bytes.toBytes(columns[11]));
                put.addColumn(CF, ATTR12_ORIGIN_AIRPORT_SEQ_ID, Bytes.toBytes(columns[12]));
                put.addColumn(CF, ATTR13_ORIGIN_CITY_MARKET_ID, Bytes.toBytes(columns[13]));
                put.addColumn(CF, ATTR14_DEST_AIRPORT_ID, Bytes.toBytes(columns[14]));
                put.addColumn(CF, ATTR15_WHEELS_ON, Bytes.toBytes(columns[15]));
                put.addColumn(CF, ATTR16_ARR_TIME, Bytes.toBytes(columns[16]));
                put.addColumn(CF, ATTR17_ARR_DELAY, Bytes.toBytes(columns[17]));
                put.addColumn(CF, ATTR18_ARR_DELAY_NEW, Bytes.toBytes(columns[18]));
                put.addColumn(CF, ATTR19_CANCELLED, Bytes.toBytes(columns[19]));
                put.addColumn(CF, ATTR20_CANCELLATION_CODE, Bytes.toBytes(columns[20]));
                put.addColumn(CF, ATTR21_AIR_TIME, Bytes.toBytes(columns[21]));
                put.addColumn(CF, ATTR22_DISTANCE, Bytes.toBytes(columns[22]));

                table.put(put);
                row_num++;
            }
        }

        reader.close();

        table.close();

        connection.close();


    }

}

