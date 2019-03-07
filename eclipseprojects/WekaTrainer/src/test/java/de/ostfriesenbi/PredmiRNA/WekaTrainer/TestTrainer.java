package de.ostfriesenbi.PredmiRNA.WekaTrainer;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URISyntaxException;

import org.junit.Test;

import picocli.CommandLine;

public class TestTrainer {

	@Test
	public void test() throws URISyntaxException {
		CommandLine.call(new App(), "-i", new File(TestTrainer.class.getClassLoader().getResource("insaexample.arff").toURI()).toString(),"weka.classifiers.trees.J48");
		CommandLine.call(new App(), "-i", new File(TestTrainer.class.getClassLoader().getResource("insaexample.arff").toURI()).toString(),"--attreval","weka.attributeSelection.CfsSubsetEval","--attrsearch","weka.attributeSelection.GreedyStepwise","weka.classifiers.trees.J48");
		CommandLine.call(new App(), "-i", new File(TestTrainer.class.getClassLoader().getResource("insaexample.arff").toURI()).toString(),"weka.classifiers.trees.RandomForest");
	}

}
