import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.TableName;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.filter.Filter;

import java.io.IOException;

public class FilterTable {
    //col family
    private static final byte[] CF = "data".getBytes();
    //cols
    private static final byte[] ATTR18_ARR_DELAY_NEW = "arr_delay_new".getBytes();
    private static final byte[] ATTR19_CANCELLED = "cancelled".getBytes();

    public static void main(String[] args) throws IOException {

        Configuration config = HBaseConfiguration.create();
        config.set("hbase.zookeper.quorum", "localhost");

        Connection connection = ConnectionFactory.createConnection(config);
        Table table = connection.getTable(TableName.valueOf("flights"));

        String custom_delay = "10.0";
        byte[] start_row = "2015-01-10".getBytes();
        byte[] stop_row = "2015-01-20".getBytes();

        Filter custom_filter = new CustomFlightFilter(Float.parseFloat(custom_delay));
        //Filter custom_filter = new CustomFlightFilter();


        Scan scan = new Scan();
        scan.addFamily(CF);
        scan.setFilter(custom_filter);
        scan.setStartRow(start_row);
        scan.setStopRow(stop_row);

        ResultScanner scanner = table.getScanner(scan);

        //make sure u have restarted HBase after any changes in CustomFilter

        for (Result r : scanner) {
            String s1 = new String(r.getRow(), "UTF-8");
            String s2 = new String(r.getValue(CF, ATTR19_CANCELLED), "UTF-8");
            String s3 = new String(r.getValue(CF, ATTR18_ARR_DELAY_NEW), "UTF-8");
            System.out.println("Key:" + s1 + " Cancelled:" + s2 + ", Delay:" + s3);
        }
        System.out.println("Done");

        scanner.close();
        table.close();
        connection.close();

        //hbase> create "flights", "data"
        //mvn package
        //mvn compile exec:java -Dexec.mainClass="FillTable"
        //../hbase-1.1.2/bin/stop-hbase.sh
        //../hbase-1.1.2/bin/start-hbase.sh
        //mvn compile exec:java -Dexec.mainClass="FilterTable"
    }
}






FILTER:
import org.apache.hadoop.hbase.Cell;
import org.apache.hadoop.hbase.exceptions.DeserializationException;
import org.apache.hadoop.hbase.filter.Filter;
import org.apache.hadoop.hbase.filter.FilterBase;
import org.apache.hadoop.hbase.filter.PageFilter;
import org.apache.hadoop.hbase.util.Bytes;

import java.io.IOException;

public class CustomFlightFilter extends FilterBase {
    private float delay;
    private boolean remove_it;

    public CustomFlightFilter() {
        super();
        this.remove_it = true;
    }

    public CustomFlightFilter(float delay) {
        this.delay = delay;
        this.remove_it = true;
    }

    @Override
    public ReturnCode filterKeyValue(Cell cell) throws IOException {
        String column = new String(
                cell.getQualifierArray(),
                cell.getQualifierOffset(),
                cell.getQualifierLength());

        String str_value = new String(
                cell.getValueArray(),
                cell.getValueOffset(),
                cell.getValueLength());

        //string № 4765 has cancelled = 0.0, and arr_delay_new is empty

        String key = new String(
                cell.getRowArray(),
                cell.getRowOffset(),
                cell.getRowLength());

        float value;

        if (delay > 0.f) {
            if (column.equals("arr_delay_new")) {

                /*if (str_value.isEmpty()) {
                    value = (float)0;
                } else {
                    value = Float.parseFloat(str_value);
                }*/

                System.out.println("key:" + key + "col:" + column + " val:" + str_value);

                if (!str_value.isEmpty() && !str_value.equals("") &&  Float.parseFloat(str_value) > 10.f) {
                    remove_it = false;
                }
                System.out.println("remove_it="+remove_it);
            }
        } else {

            if (column.equals("cancelled")) {

                value = Float.parseFloat(str_value);

                if ((float)0 != value) {
                    remove_it = false;

                }
                System.out.println("key:" + key +"remove_it="+remove_it);

            }
        }

        return ReturnCode.INCLUDE;
        //insert this cell into processed row like <cell> + <cell> + <cell> ...
        //all the cells should remain in processed row
    }

    public static Filter parseFrom(byte[] pbBytes) throws DeserializationException{
        return new CustomFlightFilter(Bytes.toFloat(pbBytes));
    }

    @Override
    public byte[] toByteArray() throws IOException {
        return Bytes.toBytes(delay);
    }

    @Override
    public boolean filterRow() throws IOException {
        return remove_it;
    }

    @Override
    public void reset() throws IOException {
        remove_it = true;
    }
}

