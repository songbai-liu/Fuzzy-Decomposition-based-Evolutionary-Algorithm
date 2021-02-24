package jmetal.operators.learning;

import java.util.HashMap;

import jmetal.util.Configuration;
import jmetal.util.JMException;

public class LearningFactory {
	
	public static Learning getLearningOperator(String name, HashMap parameters) throws JMException{
		if(name.equalsIgnoreCase("BasicLearning")){
			return new BasicLearning(parameters);
		}else if(name.equalsIgnoreCase("IntegratedLearning")){
			return new IntegratedLearning(parameters);
		}else if(name.equalsIgnoreCase("DeepLearning")){
			return new DeepLearning(parameters);
		}else{
			Configuration.logger_.severe("LearningFactory.getLearningOperator. " +
			          "Operator '" + name + "' not found ");
			      throw new JMException("Exception in " + name + ".getLearningOperator()") ;
		}
	}

}
