package de.ostfriesenbi.PredmiRNA.WekaClassify;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;
import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SerializationHelper;
import weka.core.converters.ConverterUtils.DataSource;

@Command(description="Applies the given serialized classifier to the given instances",name="WekaClassifier",mixinStandardHelpOptions=true,version="WekaClassifier 0.1")
public class App implements Callable<Void>{

	@Parameters(index="0", description="The identifier of the classifier (extends AbstractClassifier), for example: weka.classifiers.trees.J48")
	private File serializedclass;

	@Option(names={"-i","--input"},description="The .arff to read the instances from", required=true )
	private File arfffile;
	
	@Option(names={"-o","--output"},description="The .csv to store the results", required=true )
	private File csvfile;
	
	@Option(names={"-n","--name"},description="The name of the identifier field in the .arff file")
	private String name="comment";
	
	@Option(names={"-t","--threads"},description="The number of threads to use")
	private int threads=2;






	public static void main( String[] args ){
		CommandLine.call(new App(), args);
	}

	@Override
	public Void call() throws Exception {
		final DataSource ds = new DataSource(arfffile.getAbsolutePath());
		final Instances instances = ds.getDataSet();
		final ExecutorService es= Executors.newFixedThreadPool(threads);
		final Classifier c=(Classifier) SerializationHelper.read(serializedclass.getAbsolutePath());
		final Map<String, Future<Double>> results=new HashMap<>();
		
		Attribute att = instances.attribute(name);
		
		for (Instance instance : instances) {
			results.put(instance.stringValue(att), es.submit(() -> c.classifyInstance(instance)));
		}
		
		try(PrintWriter pw=new PrintWriter(csvfile)){
			pw.println(name+",realmiRNA");
			for (Entry<String, Future<Double>> result : results.entrySet()) {
				pw.println(result.getKey()+","+result.getValue().get().toString());
			}
		}
		return null;
	}
}
