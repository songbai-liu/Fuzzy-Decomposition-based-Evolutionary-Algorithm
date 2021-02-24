package jmetal.operators.mutation;

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

public class M2MMutation extends Mutation{
	private Double mutationProb = null;
	
	/**
	 * Valid solution types to apply this operator
	 */
	private static final List VALID_TYPES = Arrays.asList(
			RealSolutionType.class, ArrayRealSolutionType.class);

	public M2MMutation(HashMap<String, Object> parameters) {
		super(parameters);
		if (parameters.get("probability") != null)
			mutationProb = (Double) parameters.get("probability");
	}
	
	/**
	 * Perform the mutation operation
	 * 
	 * @param probability
	 *            Mutation probability
	 * @param solution
	 *            The solution to mutate
	 * @throws JMException
	 */
	public void doMutation(double probability, double delta, Solution solution)
			throws JMException {
		double rnd;
		double rand1, rand2,rand3;
		double y, yl, yu, val, xy;
		XReal x = new XReal(solution);
		for (int var = 0; var < solution.numberOfVariables(); var++) {
			if (PseudoRandom.randDouble() <= probability) {
				y = x.getValue(var);
				yl = x.getLowerBound(var);
				yu = x.getUpperBound(var);
				rand1 = PseudoRandom.randDouble();
				//rand2 = PseudoRandom.randDouble();
				rnd = 0.5*(1.0-Math.pow(rand1, -delta))*(rand1-0.5);
				y = y+rnd*(yu-yl);
				rand3 = PseudoRandom.randDouble();
				if(y < yl){
					y = (yl + 0.5*rand3*(x.getValue(var) - yl));
				}
				if(y > yu){
					y = (yu - 0.5*rand3*(yu -  x.getValue(var)));
				}
				x.setValue(var, y);
			}
		}
	}

	public Object execute(Object object, double delta) throws JMException {
		Solution solution = (Solution) object;

		if (!VALID_TYPES.contains(solution.getType().getClass())) {
			Configuration.logger_
					.severe("M2MMutation.execute: the solution "
							+ "type " + solution.getType()
							+ " is not allowed with this operator");

			Class cls = java.lang.String.class;
			String name = cls.getName();
			throw new JMException("Exception in " + name + ".execute()");
		} // if

		doMutation(mutationProb, delta, solution);
		return solution;
	}
	
	public Object execute(Object object) throws JMException {
		// TODO Auto-generated method stub
		return null;
	}

}
