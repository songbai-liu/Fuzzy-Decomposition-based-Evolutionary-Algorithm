package jmetal.operators.crossover;

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

public class M2MCrossover extends Crossover{
	//private static final double DEFUALT_DELTA = 0.5;
	//private double delta;
	
	/**
	 * Valid solution types to apply this operator
	 */
	private static final List VALID_TYPES = Arrays.asList(
			RealSolutionType.class, ArrayRealSolutionType.class);

	/**
	 * Constructor
	 */
	public M2MCrossover(HashMap<String, Object> parameters) {
		super(parameters);
		//delta = DEFUALT_DELTA;
	}
	
	/**
	 * Perform the crossover operation.
	 * 
	 * @param probability
	 *            Crossover probability
	 * @param parent1
	 *            The first parent
	 * @param parent2
	 *            The second parent
	 * @return An array containing the two offsprings
	 */
	public Solution doCrossover(double delta, Solution parent1,
			Solution parent2) throws JMException {

		Solution offSpring = new Solution(parent1);
		int i;
		double rand1, rand2, rand3;
		double yL, yu;
		double c;
		double rnd;
		double valueX1, valueX2;
		XReal x1 = new XReal(parent1);
		XReal x2 = new XReal(parent2);
		XReal offs = new XReal(offSpring);

		int numberOfVariables = x1.getNumberOfDecisionVariables();
		for (i = 0; i < numberOfVariables; i++) {
			valueX1 = x1.getValue(i);
			valueX2 = x2.getValue(i);
			
			yL = x1.getLowerBound(i);
			yu = x1.getUpperBound(i);
			rand1 = PseudoRandom.randDouble();
			//rand2 = PseudoRandom.randDouble();
			rnd = (1.0-Math.pow(rand1, -delta))*2.0*(rand1-0.5);
			c = valueX1 + rnd*(valueX2-valueX1);
			rand3 = PseudoRandom.randDouble();
			if(c < yL){
				c = (yL + 0.5*rand3*(offs.getValue(i) - yL));
			}
			if(c > yu){
				c = (yu - 0.5*rand3*(yu - offs.getValue(i)));
			}
			offs.setValue(i, c);
		} // if
			
		return offSpring;
	} // doCrossover

	/**
	 * Executes the operation
	 * 
	 * @param object
	 *            An object containing an array of two parents
	 * @return An object containing the offSprings
	 */
	public Object execute(Object object, double delta) throws JMException {
		Solution[] parents = (Solution[]) object;

		if (parents.length != 2) {
			Configuration.logger_
					.severe("M2MCrossover.execute: operator needs two "
							+ "parents");
			Class cls = java.lang.String.class;
			String name = cls.getName();
			throw new JMException("Exception in " + name + ".execute()");
		} // if

		if (!(VALID_TYPES.contains(parents[0].getType().getClass()) && VALID_TYPES
				.contains(parents[1].getType().getClass()))) {
			Configuration.logger_.severe("M2MCrossover.execute: the solutions "
					+ "type " + parents[0].getType()
					+ " is not allowed with this operator");

			Class cls = java.lang.String.class;
			String name = cls.getName();
			throw new JMException("Exception in " + name + ".execute()");
		} // if

		Solution offSpring;
		offSpring = doCrossover(delta, parents[0], parents[1]);

		// for (int i = 0; i < offSpring.length; i++)
		// {
		// offSpring[i].setCrowdingDistance(0.0);
		// offSpring[i].setRank(0);
		// }
		return offSpring;
	} // execute

	@Override
	public Object execute(Object object) throws JMException {
		// TODO Auto-generated method stub
		return null;
	}

}
