package jmetal.operators.learning;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import jmetal.core.Solution;
import jmetal.encodings.solutionType.ArrayRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.wrapper.XReal;

public class BasicLearning extends Learning{
	private double learningProbability_ = 0.5;
	private double learningRate_ = 0.5;
	
	private static final List VALID_TYPES = Arrays.asList(RealSolutionType.class,ArrayRealSolutionType.class);
	
	public BasicLearning(HashMap<String, Object> parameters){
		super(parameters);
		if (parameters.get("learningProbability") != null)
			learningProbability_ = (Double) parameters.get("learningProbability");
		if (parameters.get("learningRate") != null)
			learningRate_ = (Double) parameters.get("learningRate");
	}
	
	public void doLearning(double probability, 
			               Solution solution1, 
			               Solution solution2) throws JMException{
		Solution child = new Solution(solution1);
		Solution parent = new Solution(solution2);
		
		double r1,c1;
		double y1,y2,yL,yU;
		double value;
		XReal Xchild = new XReal(solution1);
		XReal Xparent = new XReal(solution2);
		XReal offspring = new XReal(child);
		
		int numberOfVariables = Xchild.getNumberOfDecisionVariables();
		
		if(PseudoRandom.randDouble() <= learningProbability_){
			for(int i=0; i<numberOfVariables; i++){
				if(PseudoRandom.randDouble() <= learningRate_){
					r1 = PseudoRandom.randDouble(0, 1);
					c1 = 1.0;
					yL = Xchild.getLowerBound(i);
					yU = Xchild.getUpperBound(i);
					value = r1*c1*(Xparent.getValue(i) - Xchild.getValue(i)) ;
					if(value < yL){
						value = yL;
					}
					if(value > yU){
						value = yU;
					}
					offspring.setValue(i, value + offspring.getValue(i));
				}//if
			}//for
		}//if
	}

	@Override
	public Object execute(Object object) throws JMException {
		Solution[] solutions = (Solution[])object;
		if(solutions.length != 2){
			Configuration.logger_.severe("BasicLearning need two solution!");
			Class cls = java.lang.String.class;
			String name = cls.getName();
			throw new JMException("Exception in" + name + ".execute()");
		}//if
		if(!VALID_TYPES.contains(solutions[0].getType().getClass()) &&
				VALID_TYPES.contains(solutions[1].getType().getClass())){
			Configuration.logger_.severe("BasicLearning.execute: the solutions" + "type" + solutions[0].getType() + 
					"is not allowed with this operator");
			
			Class cls = java.lang.String.class;
			String name = cls.getName();
			throw new JMException("Exception in" + name + ".execute()");
		}//if
		doLearning(learningProbability_,solutions[0],solutions[1]);
		return solutions;
	}

	@Override
	public Object execute(Object object, double delta) throws JMException {
		// TODO Auto-generated method stub
		return null;
	}

}
