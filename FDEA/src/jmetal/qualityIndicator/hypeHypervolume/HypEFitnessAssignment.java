package jmetal.qualityIndicator.hypeHypervolume;



import jmetal.core.Solution;
import jmetal.core.SolutionSet;

public class HypEFitnessAssignment {
	
	public void setHypEFitness(SolutionSet population, Solution reference, int k, int nrOfSamples) {
		if (reference.numberOfObjectives() <= 3)
			setExactHypEFitness(population, reference, k);
		else
			setEstimatedHypEFitness(population, reference, k, nrOfSamples);
	}
	
	
	void setExactHypEFitness(SolutionSet population, Solution reference, int k) {
		HypEFitness hy = new HypEFitness();

		int objs = reference.numberOfObjectives();

		int size = population.size();

		double[][] points = new double[size + 1][objs + 1];
		
		for (int i = 1; i <= size; i++) {
			for (int j = 1; j <= objs; j++)
				points[i][j] = population.get(i - 1).getObjective(j - 1);
		}

		double[] bounds = new double[objs + 1];
		for (int i = 1; i <= objs; i++)
			bounds[i] = reference.getObjective(i-1);


		double[] result = hy.hypeIndicatorExact(points, bounds, k);

		for (int i = 1; i <= size; i++) {
			population.get(i - 1).setFitness(result[i]);
		}

	}
	
	
	
	void setEstimatedHypEFitness(SolutionSet population,
			Solution reference, int k, int nrOfSamples) {

		HypEFitness hy = new HypEFitness();

		int objs = reference.numberOfObjectives();

		int size = population.size();

		double[][] points = new double[size][objs];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < objs; j++)
				points[i][j] = population.get(i).getObjective(j);
		}

		double lowerbound = 0;
		double upperbound = reference.getObjective(0);

		double[] result = hy.hypeIndicatorSampled(points, lowerbound, upperbound, k,
				nrOfSamples);

		for (int i = 0; i < size; i++) {
			population.get(i).setFitness(result[i]);
		}

	}
	

}
