package jmetal.metaheuristics.moeaddu;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PairData;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Comparator;
import java.util.StringTokenizer;
import java.util.Vector;


// This is the implementation of MOEA/D-DU without normalization


public class MOEADDU extends Algorithm {

	private int populationSize_;  // population size

	private SolutionSet population_;  // current population

	double[] z_;   // ideal point

	double[][] lambda_;   // well-distributed weight vectors

	int T_;   // neighborhood size
	int K_;   // parameter K 

	int[][] neighborhood_;  // indexes of neighborhood

	double delta_;   // probability to choose the mating solution from the neighborhood

	
	String functionType_;  // type of aggregation function
	
	int generations_;   // objective function evaluations
	
	
	/**
	 * Operators
	 */
	Operator crossover_;  // crossover operator
	Operator mutation_;   // mutation operator

	int div1_;
	int div2_;
	
	

	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public MOEADDU(Problem problem) {
		super(problem);

		functionType_ = "_TCHE";

	} // DMOEA

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxGenerations;

		
		generations_ = 0;
		
		maxGenerations = ((Integer) this.getInputParameter("maxGenerations"))
				.intValue();  
		
		div1_ = ((Integer) this.getInputParameter("div1")).intValue();
		div2_ = ((Integer) this.getInputParameter("div2")).intValue();
		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,
				problem_.getNumberOfObjectives());
		lambda_ = vg.getVectors();

		populationSize_ = vg.getVectors().length;

		population_ = new SolutionSet(populationSize_);

		T_ = ((Integer) this.getInputParameter("T")).intValue();
		K_ = ((Integer) this.getInputParameter("K")).intValue();
		delta_ = ((Double) this.getInputParameter("delta")).doubleValue();

		neighborhood_ = new int[populationSize_][T_];

		z_ = new double[problem_.getNumberOfObjectives()];

		crossover_ = operators_.get("crossover"); // SBX crossover
		mutation_ = operators_.get("mutation"); // polynomial mutation


		initNeighborhood(); // initialize neighborhood

		
		initPopulation();  // initialize population


		initIdealPoint();  // initialize ideal point

		// STEP 2. Update
		do {
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);

			for (int i = 0; i < populationSize_; i++) {
				int n = permutation[i]; // or int n = i;
	
				int type;
				double rnd = PseudoRandom.randDouble();

				// STEP 2.1. Mating selection based on probability
				if (rnd < delta_) // if (rnd < realb)
				{
					type = 1; // neighborhood
				} else {
					type = 2; // whole population
				}
				
				
				Vector<Integer> p = new Vector<Integer>();
				matingSelection(p, n, 1, type);

				// STEP 2.2. Reproduction
				Solution child;
				Solution[] parents = new Solution[2];

				parents[0] = population_.get(p.get(0));
				parents[1] = population_.get(n);


				Solution[] offSpring = (Solution[]) crossover_.execute(parents);

				child = offSpring[0];

				// Apply mutation
				mutation_.execute(child);

				// Evaluation
				problem_.evaluate(child);
					
				updateIdealPoint(child);

				updateProblem(child);
			} // for
		//	one generation
			generations_++;
		} while (generations_ < maxGenerations);

		
		NondominatedRanking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);
	}



	void updateProblem(Solution indiv) {
		PairData[] pairs = new PairData[populationSize_];

		for (int i = 0; i < populationSize_; i++) {

			double dist = getPerpendicularDistance(indiv, lambda_[i]);
			pairs[i] = new PairData(i, dist);
		}

		Arrays.sort(pairs);
		
		
		for (int i = 0; i < K_; i++) {
			
			int n = i;
			int k = pairs[n].getID();

			double f1, f2;

			f1 = fitnessFunction(population_.get(k), lambda_[k]);
			f2 = fitnessFunction(indiv, lambda_[k]);

			if (f2 < f1) {
				population_.replace(k, new Solution(indiv));
				return;  
			}

		}

	} // updateProblem
	
	

	double getPerpendicularDistance(Solution individual, double[] lambda) {
		double d1, d2, nl;

		d1 = d2 = nl = 0.0;

		for (int i = 0; i < individual.getNumberOfObjectives(); i++) {
			d1 += (individual.getObjective(i) - z_[i]) * lambda[i];

			nl += (lambda[i] * lambda[i]);
		}
		nl = Math.sqrt(nl);
		d1 = Math.abs(d1) / nl;

		for (int i = 0; i < individual.getNumberOfObjectives(); i++) {

			d2 += ((individual.getObjective(i) - z_[i]) - d1 * (lambda[i] / nl))
					* ((individual.getObjective(i) - z_[i]) - d1
							* (lambda[i] / nl));
		}
		d2 = Math.sqrt(d2);

		return d2;

	}
	
	public void initNeighborhood() {
		double[] x = new double[populationSize_];
		int[] idx = new int[populationSize_];

		for (int i = 0; i < populationSize_; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize_; j++) {
				x[j] = Utils.distVector(lambda_[i], lambda_[j]);
				idx[j] = j;
			} // for

			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize_, T_);
			// minfastsort(x,idx,population.size(),niche);

			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
		} // for
	} // initNeighborhood

	/**
* 
*/
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);

			//generations_++;
			population_.add(newSolution);
		} // for
	} // initPopulation

	/**
* 
*/
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			z_[i] = 1.0e+30;
		
		for (int i = 0; i < populationSize_; i++)
			updateIdealPoint(population_.get(i));
	} // initIdealPoint

	/**
* 
*/
	public void matingSelection(Vector<Integer> list, int cid, int size,
			int type) {
		// list : the set of the indexes of selected mating parents
		// cid : the id of current subproblem
		// size : the number of selected mating parents
		// type : 1 - neighborhood; otherwise - whole population
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
				// p = population[cid].table[r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection

	/**
	 * 
	 * @param individual
	 */
	void updateIdealPoint(Solution individual) {
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < z_[n])
				z_[n] = individual.getObjective(n);
		}
	} // updateReference




	
	double fitnessFunction(Solution individual, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TCHE")) {
			double maxFun = -1.0e+30;

			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				double diff = Math.abs(individual.getObjective(n) - z_[n]);

				double feval;
				if (lambda[n] == 0) {
					feval = diff / 0.000001;
				} else {
					feval = diff / lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for

			fitness = maxFun;
		} // if
		else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation // MOEAD

} // MOEAD

