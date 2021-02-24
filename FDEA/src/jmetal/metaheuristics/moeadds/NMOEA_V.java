
package jmetal.metaheuristics.moeadds;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.ranking.NondominatedRanking;

/**
 * 
 */

public class NMOEA_V extends Algorithm {
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet union_;
	
	private int populationSize_;
	int generations_;
	int maxGenerations_;
	
	Operator selection_;
	Operator crossover_;
	Operator mutation_;
	
	double[][] lamda_;
	private int div1_;
	private int div2_;

	
	private double[] zideal_; //ideal point
	private double[] znadir_;//Nadir point
	
	/**
	 * T: neighbour size
	 */
	int T_;
	/**
	 * Neighborhood
	 */
	int[][] neighborhood_;
	/**
	 * delta: probability that parent solutions are selected from neighbourhood
	 */
	double delta_;
	
	String dataDirectory_;
	
	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public NMOEA_V(Problem problem) {
		super(problem);
		zideal_ = new double[problem.getNumberOfObjectives()];
		znadir_ = new double[problem.getNumberOfObjectives()];
	}

	/**
	 * Runs the MOEA/C algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		/*
		 * step1: Basic Setting of this Algorithm
		 */
		baiscSetting();
		/*
		 * step2: Initialize the Population
		 */
		initPopulation();
		
		initUniformWeight();
		/*
		 * Enter the main loop，into the process of evolution
		 */
		while (generations_ < maxGenerations_) {
			/*
			 * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
			 */
			//generateOffspringPopulation();
			generateOffspringBasedNeighbor_DE();
			//generateOffspringBasedNeighbor();
			/*
			 * step5:Environmental Selection
			 */
		    environmentalSelection();
		    
		    generations_++;
		}
		
		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);
	} // execute
	
	/*
	 * step1: Basic Setting of this Algorithm
	 */
	public void baiscSetting(){
		generations_ = 0;
		maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations")).intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
		mutation_  = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");
		
		T_ = ((Integer) this.getInputParameter("T")).intValue();
		//T_ = (int)0.1*populationSize_;
		delta_ = ((Double) this.getInputParameter("delta")).doubleValue();
		neighborhood_ = new int[populationSize_][T_];
		lamda_ = new double[populationSize_][problem_.getNumberOfObjectives()];
		dataDirectory_  = this.getInputParameter("dataDirectory").toString();
	}
	
	/*
	 * step2: Initialize the Population
	 */
	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population_.add(newSolution);
			//generations_++;
		} // for
		estimateIdealPoint(population_);
		estimateNadirPoint(population_);
		normalizationObjective(population_);
		computeDistance(population_);
	} // initPopulation
	
	/*
	 * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
	 */
	public void generateOffspringPopulation() throws JMException{
		offspringPopulation_ = new SolutionSet(populationSize_);
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize_); i++) {
			
			// obtain parents
			/*parents = (Solution[]) selection_
					.execute(population_);*/
			parents[0] = (Solution) selection_.execute(population_);
			parents[1] = (Solution) selection_.execute(population_);
			Solution offSpring = (Solution) crossover_
					.execute(parents);
			mutation_.execute(offSpring);
			problem_.evaluate(offSpring);
			problem_.evaluateConstraints(offSpring);
			offspringPopulation_.add(offSpring);
			//generations_++;
		} // for
	}
	
	public void generateOffspringBasedNeighbor() throws JMException{
		getNeighborhood();
		offspringPopulation_ = new SolutionSet(populationSize_);
		Solution[] parents = new Solution[2];
		int[] permutation = new int[populationSize_];
		Utils.randomPermutation(permutation, populationSize_);
		for (int i = 0; i < populationSize_; i++) {
			double delta = Math.pow(1.0-(double)generations_/maxGenerations_, 0.7);
			int n = permutation[i]; // or int n = i;
			// int n = i ; // or int n = i;
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
			matingSelection(p, n, 2, type);
			parents[0] = population_.get(p.get(0));
			parents[1] = population_.get(p.get(1));
			Solution[] offSpring = (Solution[]) crossover_
					.execute(parents);
			mutation_.execute(offSpring[0]);
			problem_.evaluate(offSpring[0]);
			problem_.evaluateConstraints(offSpring[0]);
			offspringPopulation_.add(offSpring[0]);
		}
	}
	
	public void generateOffspringBasedNeighbor_DE() throws JMException{
		getNeighborhood();
		offspringPopulation_ = new SolutionSet(populationSize_);
		Solution[] parents = new Solution[3];
		int[] permutation = new int[populationSize_];
		Utils.randomPermutation(permutation, populationSize_);
		for (int i = 0; i < populationSize_; i++) {
			int n = permutation[i]; // or int n = i;
			// int n = i ; // or int n = i;
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
			matingSelection(p, n, 2, type);
			parents[0] = population_.get(p.get(0));
			parents[1] = population_.get(p.get(1));
			parents[2] = population_.get(n);
			Solution offSpring = (Solution) crossover_.execute(new Object[] {
					population_.get(n), parents });
			mutation_.execute(offSpring);
			problem_.evaluate(offSpring);
			problem_.evaluateConstraints(offSpring);
			offspringPopulation_.add(offSpring);
		}
	}

	/*
	 * step5:Environmental Selection
	 */
	public void environmentalSelection(){
		/*
		 * step5.1:Combine the Population and the Offspring Population
		 */
		union_ = ((SolutionSet) population_).union(offspringPopulation_);
		/*
		 * step5.2:Normalization the Combined Population
		 */
		//Ranking nodominatedRanking = new NondominatedRanking(union_);
		SolutionSet[] st = getStSolutionSet(union_,populationSize_);
		estimateIdealPoint(st[1]);
		estimateNadirPoint(st[1]);
		//estimateIdealPoint(union_);
		//estimateNadirPoint(union_);
		normalizationObjective(union_);
		/*
		 * step5.3:Compute the Convergence Distance of each Solution
		 */
		computeDistance(union_);
		
		SolutionSet[] clusterSets = associate(union_);
		

		
		population_.clear();
      

		bestSolutionSelection(clusterSets);	
			
	}
	
	public SolutionSet[] associate(SolutionSet solutionSet){
		SolutionSet[] clusters = new SolutionSet[populationSize_];
		int[] numb = new int[populationSize_];
		double[] mAngles = new double[populationSize_];
		for(int i=0;i<populationSize_;i++){
			clusters[i] = new SolutionSet();
			numb[i] = 0;
			mAngles[i] = Double.MIN_VALUE;
		}
		for (int k = 0; k < solutionSet.size(); k++) {

			Solution sol = solutionSet.get(k);

			//double min = calVDistance(sol, lamda_[0]);
			double min = calNormlizedVDistance(sol, lamda_[0]);
			int index = 0;

			for (int j = 1; j < lamda_.length; j++) {
				//double dist = calVDistance(sol, lamda_[j]);
				double dist = calNormlizedVDistance(sol, lamda_[j]);
				if (dist < min) {
					min = dist;
					index = j;
				}
			}
			sol.setClusterID(index);
			clusters[index].add(sol);
			numb[index]++;
		}
		int size = 0;
		int[] emptyID = new int[populationSize_];
		for(int i=0;i<populationSize_;i++){
			if(numb[i] == 0){
				emptyID[size] = i;
				size++;
			}else{
				for(int j=0;j<clusters[i].size();j++){
					for(int k=0;k<clusters[i].size();k++){
						if(k!=j){
							double ang = computeAngle(clusters[i].get(j),clusters[i].get(k));
							if(mAngles[i] < ang){
								mAngles[i] = ang;
							}
						}
					}
					
				}
			}
			
		}
       int t = size;
		while(size > 0){
			/*int id = -1;
			if(size > t/2){
				id = findMaxAngleID(mAngles);
			}else{
				id = findMaxSizeID(numb);
			}*/
			int id = findMaxSizeID(numb);
			//int id = findMaxAngleID(mAngles);
			SolutionSet[] sets = furtherDividCluster(clusters[id]);
			double angle0 = Double.MIN_VALUE;
			double angle1 = Double.MIN_VALUE;
			/*for(int i=0;i<sets[0].size();i++){
				for(int j=0;j<sets[0].size();j++){
					if(i!=j){
						double anl = computeAngle(sets[0].get(i),sets[0].get(j));
						if(angle0 < anl){
							angle0 = anl;
						}
					}
				}
			}
			
			for(int i=0;i<sets[1].size();i++){
				for(int j=0;j<sets[1].size();j++){
					if(i!=j){
						double anl = computeAngle(sets[1].get(i),sets[1].get(j));
						if(angle1 < anl){
							angle1 = anl;
						}
					}
				}
			}*/
			clusters[id] = sets[0];
			numb[id] = clusters[id].size();
			mAngles[id] = angle0;
			clusters[emptyID[size-1]] = sets[1];
			numb[emptyID[size-1]] = sets[1].size();
			mAngles[emptyID[size-1]] = angle1;
			size--;
		}
		return clusters;

	}
	
	int findMaxSizeID(int[] ind){
		int id = -1;
		int max = Integer.MIN_VALUE;
		for(int i=0;i<ind.length;i++){
			if(max < ind[i]){
				max=ind[i];
				id = i;
			}
		}
		return id;
	}
	
	int findMaxAngleID(double[] angs){
		int id = -1;
		double maxAngle = Double.MIN_VALUE;
		for(int i=0;i<angs.length;i++){
			if(maxAngle < angs[i]){
				maxAngle=angs[i];
				id = i;
			}
		}
		return id;
	}
	
	SolutionSet[] furtherDividCluster(SolutionSet solutionSet){
		List<SolutionSet> list = new <SolutionSet>ArrayList();
		for(int i=0;i<solutionSet.size();i++){
 			SolutionSet sols = new SolutionSet();
 			sols.add(solutionSet.get(i));
 			list.add(sols);
 		}
		list = new DividCluster(list).clusteringAnalysis(2);
		SolutionSet[] solSets = new SolutionSet[2];
		for(int i=0;i<2;i++){
			solSets[i] = new SolutionSet();
			for(int j=0;j<list.get(i).size();j++){
				//Solution sol = new Solution(list.get(i).get(j));
				solSets[i].add(list.get(i).get(j));
			}
		}
		return solSets;
	}
	
    public double calVDistance(Solution sol, double[] ref){
    	double d1, d2, nl;

		d1 = d2 = nl = 0.0;

		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d1 += (sol.getObjective(i) - zideal_[i]) * ref[i];
			nl += (ref[i] * ref[i]);
		}
		nl = Math.sqrt(nl);
		d1 = Math.abs(d1) / nl;
		
	
		d2 =0;
		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			
			d2 += ((sol.getObjective(i) - zideal_[i]) - d1
					* (ref[i] / nl)) * ((sol.getObjective(i) - zideal_[i]) - d1
							* (ref[i] / nl));
		}
		d2 = Math.sqrt(d2);
		

		
		return d2;
    }
    
    public double calNormlizedVDistance(Solution sol, double[] ref) {

		double ip = 0;
		double refLenSQ = 0;

		for (int j = 0; j < sol.numberOfObjectives(); j++) {

			ip += sol.getNormalizedObjective(j) * ref[j];
			refLenSQ += (ref[j] * ref[j]);
		}
		refLenSQ = Math.sqrt(refLenSQ);

		double d1 = Math.abs(ip) / refLenSQ;

		double d2 = 0;
		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d2 += (sol.getNormalizedObjective(i) - d1 * (ref[i] / refLenSQ))
					* (sol.getNormalizedObjective(i) - d1 * (ref[i] / refLenSQ));
		}
		d2 = Math.sqrt(d2);

		return d2;
	}
	
	/*
	 * Estimate the Ideal Point 
	 */
	public void estimateIdealPoint(SolutionSet solutionSet){
		for(int i=0; i<problem_.getNumberOfObjectives();i++){
			zideal_[i] = 1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) < zideal_[i]){
					zideal_[i] = solutionSet.get(j).getObjective(i);
				}
			}
			
		}
	}
	
	/*
	 * Estimate the Nadir Point 
	 */
    public void estimateNadirPoint(SolutionSet solutionSet){
    	for(int i=0; i<problem_.getNumberOfObjectives();i++){
			znadir_[i] = -1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) > znadir_[i]){
					znadir_[i] = solutionSet.get(j).getObjective(i);
				}
			}
			
		}
	}
	
    /*
     * Normalization
     */
	public void normalizationObjective(SolutionSet solutionSet){
		for(int i=0; i<solutionSet.size(); i++){
			Solution sol = solutionSet.get(i);
			
			for(int j=0; j<problem_.getNumberOfObjectives(); j++){
				double val = 0.0;
				val = (sol.getObjective(j) - zideal_[j])/(znadir_[j]-zideal_[j]);
				//val = (sol.getObjective(j) - zideal_[j]);
				sol.setNormalizedObjective(j, sol.getObjective(j));
				//sol.setNormalizedObjective(j, val);
			}
		}
	}
	
	 /*
     * Compute the Convergence Distance of each Solutions Which use the distance of 
     * each solution to the Ideal Point
     */
    public void computeDistance(SolutionSet solutionSet){
    	for(int i=0; i<solutionSet.size(); i++){
    		Solution sol = solutionSet.get(i);
    		double normDistance = 0.0;
    		double sumValue = 0.0;
    		double nadirDistance = 0.0;
    		for(int j=0; j<problem_.getNumberOfObjectives(); j++){
    			normDistance += sol.getNormalizedObjective(j) * sol.getNormalizedObjective(j);
    			sumValue +=  sol.getNormalizedObjective(j);
    			nadirDistance += (sol.getNormalizedObjective(j)-1.0)*(sol.getNormalizedObjective(j)-1.0);
    		}
    		normDistance = Math.sqrt(normDistance);
    		nadirDistance = Math.sqrt(nadirDistance);
    		
    		sol.setDistanceToIdealPoint(normDistance);
    		sol.setDistanceToNadirPoint(nadirDistance);
    		sol.setSumValue(sumValue);
    	}
    }
    
	public double weightSumValue(Solution solution){
		double value = 0.0;
		for(int i=0; i<problem_.getNumberOfObjectives();i++){
			value += solution.getNormalizedObjective(i);
		}
		return value;
	}
	public double weightSumValue(Solution r,Solution s){
		double value = 0.0;
		double[] lamda = new double[problem_.getNumberOfObjectives()];
		double sum = weightSumValue(r);
		for(int i=0; i<problem_.getNumberOfObjectives();i++){
			lamda[i] = r.getNormalizedObjective(i)/sum;
			value += lamda[i]*s.getNormalizedObjective(i);
		}
		return value;
	}
	
	public double computeChebyshev(Solution r, Solution s){
		double fitness;
		fitness = 0.0;

		double maxFun = -1.0e+30;

		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			double diff = Math.abs(s.getNormalizedObjective(n));

			double feval;
			if (r.getNormalizedObjective(n) == 0) {
				feval = diff / 0.000001;
			} else {
				feval = diff / r.getNormalizedObjective(n);
			}
			if (feval > maxFun) {
				maxFun = feval;
			}
		} // for

		fitness = maxFun;
		return fitness;
	}
	
	public double computeASFFitness(Solution r,Solution s){
		double fitness = -1.0e+30;
		int objectiveSize = r.getNumberOfObjectives();
		double[] lambda = new double[objectiveSize];
		double sumValue = 0.0;
		for(int i=0; i<objectiveSize; i++){
			sumValue += r.getNormalizedObjective(i);
		}
		for(int j=0; j<objectiveSize; j++){
			lambda[j] = r.getNormalizedObjective(j)/sumValue;
			if(lambda[j] == 0){
				lambda[j] = 0.000001;
			}
		}
		for(int k=0; k<objectiveSize; k++){
			double sb = s.getNormalizedObjective(k)/lambda[k];
			if(fitness < sb){
				fitness = sb;
			}
		}
		return fitness;
	}
	
	 /*
     * Compute the angle value between Solution1 and Solution2
     */
	public double computeAngle(Solution s1, Solution s2){
		double angle = 0.0;
		double distanceToidealPoint1 = s1.getDistanceToIdealPoint();
		double distanceToidealPoint2 = s2.getDistanceToIdealPoint();
		double innerProduc = 0.0; 
		for(int i=0; i<problem_.getNumberOfObjectives(); i++){
			innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
		}
		double value = innerProduc/(distanceToidealPoint1*distanceToidealPoint2);
		if(value > 1.0){
			value = 1.0;
		}
		angle = Math.acos(Math.abs(value));
		//System.out.println(Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2)));
		return angle;
	}//computeAngle
	
	public double computePBIFitness(Solution r,Solution s){
		double fitness = 0.0;
		double d1,d2,norm;
		double objectiveSize = r.getNumberOfObjectives();
		d1 = d2 = norm = 0.0;
		for(int i=0; i<objectiveSize; i++){
			d1 += r.getNormalizedObjective(i)*s.getNormalizedObjective(i);
			norm += r.getNormalizedObjective(i)*r.getNormalizedObjective(i);
		}
		norm = Math.sqrt(norm);
		d1 = Math.abs(d1)/norm;
		for(int j=0; j<objectiveSize; j++){
			d2 += (s.getNormalizedObjective(j)-d1*(r.getNormalizedObjective(j)/norm))
					*(s.getNormalizedObjective(j)-d1*(r.getNormalizedObjective(j)/norm));
			//d2 += s.getNormalizedObjective(j)*s.getNormalizedObjective(j);
		}
	/*	d2 = d2 - Math.pow(d1, 2.0);
		if(d2 < 0){
			//System.out.println(d2);
			d2 = 0;
		}*/
		d2 = Math.sqrt(d2);
		fitness = d1 + 3.0*d2;
		return fitness;
	}
    
    public SolutionSet[] getStSolutionSet(SolutionSet ss,int size) {
		
		SolutionSet[] sets = new SolutionSet[2];
		Ranking ranking = new NondominatedRanking(ss);

		int remain = size;
		int index = 0;
		SolutionSet front = null;
		SolutionSet mgPopulation = new SolutionSet();
		front = ranking.getSubfront(index);
		sets[0] = front;
		while ((remain > 0) && (remain >= front.size())) {

			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			} // for

			// Decrement remain
			remain = remain - front.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		}
		if (remain > 0) { // front contains individuals to insert
			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			}
		}

		sets[1] = mgPopulation;

		return sets;
	 }
 
    public void bestSolutionSelection(SolutionSet[] solutionSets){
    	for(int i=0;i<solutionSets.length;i++){
    		double minValue = Double.MAX_VALUE;
    		int ids = -1;
    		if(solutionSets[i].size() == 0){
    			System.out.println("solutionSets[i].size() == 0");
    			System.exit(0);
    		}
            Solution sol1 = solutionSets[i].getCentroidVector();
    		for(int j=0;j<solutionSets[i].size();j++){
    			//double value = computeChebyshev(sol1,solutionSets[i].get(j));
    			double value = weightSumValue(sol1,solutionSets[i].get(j));
    			//if(minValue > solutionSets[i].get(j).getSumValue()){
    			if(minValue > value){	
    				//minValue = solutionSets[i].get(j).getSumValue();
    				minValue = value;
    				ids = j;
    			}
    		}
    		population_.add(solutionSets[i].get(ids));
    	}
    }
    
    /**
     * 
     */
  	public void getNeighborhood() {
  		double[] x = new double[populationSize_];
  		int[] idx = new int[populationSize_];

  		for (int i = 0; i < populationSize_; i++) {
  			// calculate the distances based on Solutions
  			for (int j = 0; j < populationSize_; j++) {
  				x[j] = computeAngle(population_.get(i), population_.get(j));
  				idx[j] = j;
  			} // for

  			// find 'niche' nearest neighboring solutions
  			Utils.minFastSort(x, idx, populationSize_, T_);
  			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
  		} // for
  	} // initNeighborhood
  	
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

  			// if (flag) list.push_back(p);
  			if (flag) {
  				list.addElement(p);
  			}
  		}
  	} // matingSelection
  	
  	public void initUniformWeight() { // init lambda vectors
  		/*div1_ = ((Integer) this.getInputParameter("div1")).intValue();

		div2_ = ((Integer) this.getInputParameter("div2")).intValue();
		
		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,
				problem_.getNumberOfObjectives());
		lamda_ = vg.getVectors();*/
		
		String dataFileName;
		dataFileName = "W" + problem_.getNumberOfObjectives() + "D_"
				+ populationSize_ + ".dat";

		try {
			// Open the file  修改过
			FileInputStream fis = new FileInputStream(dataDirectory_ 
					+ dataFileName);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);

			int i = 0;
			int j = 0;
			String aux = br.readLine();
			while (aux != null) {
				StringTokenizer st = new StringTokenizer(aux);
				j = 0;
				while (st.hasMoreTokens()) {
					double value = (new Double(st.nextToken())).doubleValue();
					lamda_[i][j] = value;
					j++;
				}
				aux = br.readLine();
				i++;
			}
			br.close();
		} catch (Exception e) {
			
			//修改过
			System.out
					.println("initUniformWeight: failed when reading for file: "
							+ dataDirectory_ + dataFileName);
			e.printStackTrace();
		}

	} // initUniformWeight
    
} // NSGA-II
