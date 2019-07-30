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
	
        @Option(names={"-c","--classatt"},description="The name of the class attribute, if not set, the second last one will be used")
        private String classattname;




	public static void main( String[] args ){
		CommandLine.call(new App(), args);
	}

	@Override
	public Void call() throws Exception {
		final DataSource ds = new DataSource(arfffile.getAbsolutePath());
		final Instances instances = ds.getDataSet();
		final Classifier c=(Classifier) SerializationHelper.read(serializedclass.getAbsolutePath());
		
		//Attribute att = instances.attribute(name);
		
                if(instances.classIndex() ==-1){
                        instances.setClassIndex(instances.numAttributes() - 2);
                }
                if(classattname!=null){
                        final Attribute classatt=instances.attribute(classattname);
                        if(classatt==null){
                                throw new RuntimeException("Could not find the specified class attribute!");
                        }
                        instances.setClass(classatt);
                }

		
		try(PrintWriter pw=new PrintWriter(csvfile)){
			pw.println(name+","+classattname);
                        int i=0;
                        for (Instance instance : instances) {
				pw.print(i++);
				pw.print(",");
				pw.println(c.classifyInstance(instance));
			}
		}
		return null;
	}
}
