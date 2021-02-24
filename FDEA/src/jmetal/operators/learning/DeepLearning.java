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

public class DeepLearning extends Learning{
	private double learningProbability_ = 0.5;
	private double learningRate_ = 0.5;
	
	private static final List VALID_TYPES = Arrays.asList(RealSolutionType.class,ArrayRealSolutionType.class);
	
	public DeepLearning(HashMap<String, Object> parameters){
		super(parameters);
		if (parameters.get("learningProbability") != null)
			learningProbability_ = (Double) parameters.get("learningProbability");
		if (parameters.get("learningRate") != null)
			learningRate_ = (Double) parameters.get("learningRate");
	}
	
	public void doLearning(double probability, 
			               Solution solution1, 
			               Solution solution2,
			               Solution solution3,
			               Solution solution4) throws JMException{
		Solution child = new Solution(solution1);
		Solution parent = new Solution(solution2);
		Solution teacher1 = new Solution(solution3);
		Solution teacher2 = new Solution(solution4);
		
		
		double r1,c1,r2,c2,r3,c3;
		double y1,y2,yL,yU;
		double value;
		XReal Xchild = new XReal(solution1);
		XReal Xparent = new XReal(solution2);
		XReal Xteacher1 = new XReal(solution3);
		XReal Xteacher2 = new XReal(solution4);
		XReal offspring = new XReal(child);
		
		int numberOfVariables = Xchild.getNumberOfDecisionVariables();
		
		if(PseudoRandom.randDouble() <= learningProbability_){
			for(int i=0; i<numberOfVariables; i++){
				if(PseudoRandom.randDouble() <= learningRate_){
					r1 = PseudoRandom.randDouble(0, 1);
					r2 = PseudoRandom.randDouble(0, 1);
					r3 = PseudoRandom.randDouble(0, 1);
					c1 = PseudoRandom.randDouble(1.5, 2.5);
					c2 = PseudoRandom.randDouble(1.5, 2.5);
					c3 = PseudoRandom.randDouble(1.5, 2.5);
					double cons = constrictionCoefficient(c1,c2,c3);
					yU = (Xchild.getUpperBound(i) - Xchild.getLowerBound(i))/2.0;
					yL = -yU;
					value = r1*c1*(Xparent.getValue(i) - Xchild.getValue(i)) +
							r2*c2*(Xteacher1.getValue(i) - Xchild.getValue(i)) +
							r2*c2*(Xteacher2.getValue(i) - Xchild.getValue(i));
					value *= cons;
					if(value < yL){
						value = yL;
					}
					if(value > yU){
						value = yU;
					}
					double fValue = value + offspring.getValue(i);
					if(fValue < Xchild.getLowerBound(i)){
						fValue = Xchild.getLowerBound(i);
					}
					if(fValue > Xchild.getUpperBound(i)){
						fValue = Xchild.getUpperBound(i);
					}
					offspring.setValue(i, fValue);
				}//if
			}//for
		}//if
	}

	@Override
	public Object execute(Object object) throws JMException {
		Solution[] solutions = (Solution[])object;
		if(solutions.length != 4){
			Configuration.logger_.severe("IntegratedLearning need three solution!");
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
		doLearning(learningProbability_,solutions[0],solutions[1],solutions[2],solutions[3]);
		return solutions;
	}
	
	 private double constrictionCoefficient(double c1, double c2 ,double c3) {
		double rho = c1 + c2 +c3 ;
		 //rho = 1.0 ;
		if (rho <= 6) {
		  return 1.0;
		} else {
		  return 4 / (2 - rho - Math.sqrt(Math.pow(rho, 2.0) - 6.0 * rho));
//		  return -1.0;
	    }
	} // constrictionCoefficient

	@Override
	public Object execute(Object object, double delta) throws JMException {
		// TODO Auto-generated method stub
		return null;
	}
}
