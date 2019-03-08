package de.ostfriesenbi.PredmiRNA.WekaTrainer;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.Callable;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;
import weka.attributeSelection.ASEvaluation;
import weka.attributeSelection.ASSearch;
import weka.attributeSelection.CfsSubsetEval;
import weka.attributeSelection.GreedyStepwise;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.classifiers.evaluation.output.prediction.AbstractOutput;
import weka.classifiers.evaluation.output.prediction.CSV;
import weka.classifiers.evaluation.output.prediction.HTML;
import weka.classifiers.meta.AttributeSelectedClassifier;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SerializationHelper;
import weka.core.Utils;
import weka.core.converters.CSVSaver;
import weka.core.converters.ConverterUtils;
import weka.core.converters.ConverterUtils.DataSource;

@Command(description="Trains and evaluates a model using the given algorithm",name="WekaTrainer",mixinStandardHelpOptions=true,version="WekaTrainer 0.1")
public class App implements Callable<Void>{

	@Parameters(index="0", description="The identifier of the classifier (extends AbstractClassifier), for example: weka.classifiers.trees.J48")
	private String classifier;

	@Option(names={"-i","--input"},description="The .arff to read the instances from", required=true )
	private File arfffile;
	@Option(names={"-c","--classatt"},description="The name of the class attribute, if not set, the last one will be used")
	private String classattname;
	@Option(names={"--classioptions"},description="Options to pass to the classifier, seperated by ,")
	private String classoptions;
	
	@Option(names="--attreval",description="The class to evaluate the attributes, also needs the searcher (extends ASEvaluation), for example: weka.attributeSelection.CfsSubsetEval ")
	private String attributeeval;
	@Option(names="--attrevaloptions",description="Options to pass to the attribute evaluator, also needs the eval, seperated by ,")
	private String attributeevaloptions;
	
	@Option(names={"--attrsearch"},description="The class to search the attribute combinations (extends ASSearch), for example: weka.attributeSelection.GreedyStepwise")
	private String attributesearch;
	@Option(names="--attrsearchoptions",description="Options to pass to the attribute searcher, seperated by ,")
	private String attributesearchoptions;
	
	
	@Option(names={"-S","--seed"},description="The name of the class attribute")
	private long seed=1;
	@Option(names={"-f","--folds"},description="The number of folds to use in the cross fold validation process.")
	private int folds=10;
	@Option(names={"-s","--statscsvfile"},description="File to store the evaluation statistics as a CSV file")
	private File csvfile;
	@Option(names={"-o","--outputclassifier"},description="File to store the serialized classifier")
	private File serializedclass;
	@Option(names={"-t","--thresholdfile"},description="CSV file with the threshhold data")
	private File threshholdfile;





	public static void main( String[] args ){
		CommandLine.call(new App(), args);
	}

	@Override
	public Void call() throws Exception {
		final DataSource ds = new DataSource(arfffile.getAbsolutePath());
		Instances instances = ds.getDataSet();
		if(instances.classIndex() ==-1){
			instances.setClassIndex(instances.numAttributes() - 1);
		}
		if(classattname!=null){
			final Attribute classatt=instances.attribute(classattname);
			if(classatt==null){
				throw new RuntimeException("Could not find the specified class attribute!");
			}
			instances.setClass(classatt);
		}
		String[] options = null;
		if(classoptions!=null)
			options=classoptions.split(",");
		Classifier ac=AbstractClassifier.forName(classifier, options);
		final Random rand=new Random(seed);
		if(attributeeval!=null&&attributesearch!=null){
			AttributeSelectedClassifier asc = new AttributeSelectedClassifier();
			String[] evaloptions = null;
			if(attributeevaloptions!=null)
				evaloptions=attributeevaloptions.split(",");
			ASEvaluation eval = ASEvaluation.forName(attributeeval, evaloptions);
			String[] searchoptions = null;
			if(attributesearchoptions!=null)
				searchoptions=attributesearchoptions.split(",");
			ASSearch search = ASSearch.forName(attributesearch, searchoptions);
			asc.setClassifier(ac);
			asc.setEvaluator(eval);
			asc.setSearch(search);
			ac=asc;
		}
		final Evaluation eval = new Evaluation(instances);
		final AbstractOutput csv=new CSV();
		csv.setBuffer(new StringBuffer());
		eval.crossValidateModel(ac, instances, folds, rand,csv);
		eval.setMetricsToDisplay(Evaluation.getAllEvaluationMetricNames());
		final StringBuffer buffer=new StringBuffer();
		buffer.append("Cross Validation:\n");
		buffer.append(eval.toSummaryString(true));
		buffer.append(eval.toMatrixString());
		buffer.append(eval.toClassDetailsString());
		if(csvfile!=null){
			try(BufferedWriter bw= new BufferedWriter(new FileWriter(csvfile))){

			}
		}
		System.out.println(buffer);
		if(serializedclass!=null){
			ac.buildClassifier(instances);
			SerializationHelper.write(serializedclass.getAbsolutePath(), ac);
		}
		if(threshholdfile!=null){
			try(BufferedOutputStream bw= new BufferedOutputStream(new FileOutputStream(threshholdfile))){
				ThresholdCurve tc = new ThresholdCurve();
				Instances tcins=tc.getCurve(eval.predictions());
				tcins.deleteWithMissing(tcins.attribute("Lift"));;
				CSVSaver cs = new CSVSaver();
				cs.setDestination(bw);
				cs.setInstances(tcins);
				cs.writeBatch();
			}
		}
		return null;
	}
}
