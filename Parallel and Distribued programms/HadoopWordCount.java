//Параллельные и распределенные программы
//hadoop lab1
//WordCountApp.java

import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.IntWritable;

public class WordCountApp {
    public static void main(String[] args) throws Exception { 
        if (args.length != 2) { 
            System.err.println("Usage: WordCountApp <input path> <output path>"); 
            System.exit(-1); 
        } 
        Job job = Job.getInstance(); 
        job.setJarByClass(WordCountApp.class); 
        job.setJobName("Word count"); 
        FileInputFormat.addInputPath(job, new Path(args[0]));
        FileOutputFormat.setOutputPath(job, new Path(args[1]));
        job.setMapperClass(WordMapper.class); 
        job.setReducerClass(WordReducer.class); 
        job.setOutputKeyClass(Text.class); 
        job.setOutputValueClass(IntWritable.class); 
        job.setNumReduceTasks(2); 
        System.exit(job.waitForCompletion(true) ? 0 : 1); 
    } 
} 




//WordMapper.java

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;

public class WordMapper extends Mapper<LongWritable, Text, Text, IntWritable> {
    @Override
    protected void map(LongWritable key, Text value, Context context) throws IOException, InterruptedException {
        String line = value.toString();

        line = line.toLowerCase();

        String[] words = line.split("[^a-zA-Z]"); //не забыть русский текст через а-яА-Я

        for (String word : words) {
            context.write(new Text(word), new IntWritable(1));
        }
    }
}




//WordReducer.java

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import java.io.IOException;
import java.util.Iterator;

public class WordReducer extends Reducer<Text, IntWritable, Text, LongWritable> {
    @Override
    protected void reduce(Text key, Iterable<IntWritable> values, Context context) throws IOException, InterruptedException {
        long count=0;
        Iterator iter = values.iterator();
        while(iter.hasNext()) {
            iter.next();
            count++;
        }
        context.write(key, new LongWritable(count));
    }
}




//на выходе 2 txtшника в outputе и пустой success

//очистить вывод в хадупе через hadoop -fs -rmr output
//затем заново вывести через copytolocal

