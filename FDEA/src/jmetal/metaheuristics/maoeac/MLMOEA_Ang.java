package jmetal.metaheuristics.maoeac;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;
import jmetal.util.Niching;
import jmetal.util.Permutation;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.comparators.GaussianLocalDensityComparator;
import jmetal.util.comparators.LocalDensityComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.ranking.ThetaRanking;

public class MLMOEA_Ang extends Algorithm{
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet union_;
	
	private int populationSize_;
	int generations_;
	int maxGenerations_;
	
	Operator selection_;
	Operator crossover_;
	Operator mutation_;
	
	private double[] zideal_; //ideal point
	private double[] znadir_;//Nadir point
	
	double[][] lambda_; // reference points
	
	//private double dtheata;
	
	public MLMOEA_Ang(Problem problem) {
		super(problem);
		zideal_ = new double[problem.getNumberOfObjectives()];
		znadir_ = new double[problem.getNumberOfObjectives()];
	}
	
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		
		//step1: 基本设置
		baiscSetting();
		//step2: 初始化种群
		initPopulation();
		//step3:主循环，进化过程
		while (generations_ < maxGenerations_) {
			//step4:Mating Selection
			//step5:Recombination
			generateOffspringPopulation();
			////step6:Environmental Selection
		    environmentalSelection();
		    
		    generations_++;
		}
		
		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);
	}
	
	public void baiscSetting(){
		generations_ = 0;
		maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations")).intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
		mutation_  = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");
	}
	
	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population_.add(newSolution);
		} // for
	} // initPopulation
	
	public void generateOffspringPopulation() throws JMException{
		offspringPopulation_ = new SolutionSet(populationSize_);
		Solution[] parents = new Solution[2];
		for (int i = 0; i < populationSize_; i++) {
			parents = (Solution[]) selection_.execute(population_);

			Solution[] offSpring = (Solution[]) crossover_
					.execute(parents);

			mutation_.execute(offSpring[0]);
			mutation_.execute(offSpring[1]);

			problem_.evaluate(offSpring[0]);
			problem_.evaluateConstraints(offSpring[0]);
			offspringPopulation_.add(offSpring[0]);
		}
	}
	
	public void environmentalSelection(){
		union_ = ((SolutionSet) population_).union(offspringPopulation_);
		//population_.clear();
		estimateIdealPoint(union_);
		estimateNadirPoint(union_);
		normalizationObjective(union_);
		computeDistanceToIdealPoint(union_);
		
	    computeLocalDensity(union_);
	    computeCenterDistance(union_);
	    
	    lambda_ = clusteringAnalysis(union_);
	    SolutionSet[] sets = getParetoFronts();
		
		SolutionSet firstFront = sets[0];   // the first non-dominated front
		SolutionSet stPopulation = sets[1]; // the population used in theta-non-dominated ranking
	    getNextPopulation(stPopulation);
	}
	
	/*
	 * Ideal Point 的估计
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
	 * Nadir point 的估计
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
	
	public void normalizationObjective(SolutionSet solutionSet){
		for(int i=0; i<solutionSet.size(); i++){
			Solution sol = solutionSet.get(i);
			
			for(int j=0; j<problem_.getNumberOfObjectives(); j++){
				double val = 0.0;
				val = (sol.getObjective(j) - zideal_[j])/(znadir_[j]-zideal_[j]);
				//val = (sol.getObjective(j) - zideal_[j]);
				sol.setNormalizedObjective(j, val);
			}
		}
	}
	
	  /*
     * 求每个个体到理想点的距离
     */
    public void computeDistanceToIdealPoint(SolutionSet solutionSet){
    	for(int i=0; i<solutionSet.size(); i++){
    		Solution sol = solutionSet.get(i);
    		double normDistance = 0.0;
    		for(int j=0; j<problem_.getNumberOfObjectives(); j++){
    			normDistance += sol.getNormalizedObjective(j) * sol.getNormalizedObjective(j);
    		}
    		normDistance = Math.sqrt(normDistance);
    		
    		sol.setDistanceToIdealPoint(normDistance);
    	}
    }
	
	 /*
     * 求两个个体之间的角度值
     */
	public double computeAngle(Solution s1, Solution s2){
		double angle = 0.0;//所求两个向量的角度
		double distanceToidealPoint1 = s1.getDistanceToIdealPoint();//S1到理想点的距离
		double distanceToidealPoint2 = s2.getDistanceToIdealPoint();//S2到理想点的距离
		double innerProduc = 0.0; //两个向量的内积
		for(int i=0; i<problem_.getNumberOfObjectives(); i++){
			innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
		}
		angle = Math.acos(Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2)));
		return angle;
	}//computeAngle
	
	public void computeLocalDensity(SolutionSet union){
		double dc = (Math.PI*problem_.getNumberOfObjectives())/(8*populationSize_);
		double[] ld = new double[union.size()];
		
		for(int i=0;i<union.size();i++){
			ld[i] = 0;
			Solution s1 = union.get(i);
			for(int j=0;j<union.size();j++){
				Solution s2 = union.get(j);
				if(i != j){
					double angle = computeAngle(s1,s2);
					if (angle <= dc){
						ld[i] += Math.pow(Math.E, -Math.pow(angle/dc, 2));					
					}					
				}
			}
		}
		for(int s=0;s<union.size();s++){
			union.get(s).setGaussianLocalDensity(ld[s]);
		}
	}
	
	public void computeCenterDistance(SolutionSet union){
		double[] cd = new double[union.size()];
		union.sort(new GaussianLocalDensityComparator());
		for(int i=1;i<union.size();i++){
			double minAng = 1.0e+30;;
			Solution s1 = union.get(i);
			for(int j=0;j<i;j++){
				Solution s2 = union.get(j);
				double angle = computeAngle(s1,s2);
				if (angle < minAng){
					minAng = angle;					
				}					
			}
			cd[i] = minAng;
		}
		double maxAng = Double.MIN_VALUE;
		for(int k=1;k<union.size();k++){
			if(maxAng < cd[k]){
				maxAng = cd[k];
			}
		}
		cd[0] = maxAng;
		
		for(int s=0;s<union.size();s++){
			union.get(s).setCenterDistance(cd[s]);
		}
	}
	
	public double[][] clusteringAnalysis(SolutionSet union){
		double[] p = new double[union.size()];
		double[] q = new double[union.size()];
		
		SolutionSet sols = new SolutionSet(union.size());
		
		  double sump = 0.0;
		    double avep = 0.0;
		    double sumq = 0.0;
		    double aveq = 0.0;
		    for(int i=0;i<union.size();i++){
		    	p[i] = union.get(i).getGaussianLocalDensity();
		    	q[i] = union.get(i).getCenterDistance();
		    	sump += p[i];
		    	sumq += q[i];
		    }
		    avep = (double)sump/union.size();
		    aveq = sumq/union.size();
		    boolean[] isAdded = new boolean[union.size()];
		    //(p[j] ==1 && q[j] > aveq) || 
		    for(int j=0;j<union.size();j++){
		    	/*if((p[j] == 0 && q[j] > 1.0*aveq) || (p[j] > 1.5*avep && q[j] > 1.5*aveq)){
		    		sols.add(union.get(j));
		    		isAdded[j] = true;
		    	}*/
		    	union.get(j).setCrowdingDistance(p[j]*q[j]);
		    }
		   /* union.sort(new CrowdingDistanceComparator());
		    for(int s=0;s<3*populationSize_/4;s++){
		    	sols.add(union.get(s));
		    	isAdded[s] = true;
		    }*/
		    System.out.print(sols.size()+" ");
		    if(sols.size() > populationSize_){
		    	deleteSolutionBasedAngle(sols,populationSize_);
		    }else if(sols.size() < populationSize_){
		    	int remain = populationSize_ - sols.size();
		    	SolutionSet sos = new SolutionSet(union.size());
		    	for(int t=0;t<union.size();t++){
		    		if(!isAdded[t]){
		    			sos.add(union.get(t));
		    		}
		    	}
		    	addSolutionBasedAngle(sols,sos,remain);
		    }
		    double[][] lamda = new double[populationSize_][problem_.getNumberOfObjectives()];
		   for(int s=0; s<populationSize_;s++){
			  for(int m=0;m<problem_.getNumberOfObjectives();m++){
				  lamda[s][m] = sols.get(s).getNormalizedObjective(m);
			  }
		   }
		   return lamda;
	}
	
	/*
     * 不用递归的方法，用循环的方法，虽然看上去很复杂，但是时间复杂度降低很多
     * */
    public void deleteSolutionBasedAngle(SolutionSet solutionSet,int solutionSize){
    	int size = solutionSet.size();
		double minAngle = Double.MAX_VALUE;
		double angle = 0.0;
		double[] minAngles = new double[size];//每个个体与之最近个体之间的角度
		int[] minIndexs = new int[size];//每个个体与之最近的个体
		int[] index = new int[2];
		index[0] = index[1] = -1;
		for(int i=0; i<size; i++){
			for(int j=0;j<size;j++){
				if(i != j){
					angle = computeAngle(solutionSet.get(i),solutionSet.get(j));
					if(minAngle > angle){
						   minAngle = angle;
						   index[0] = i;
						   index[1] = j;
					}
				}
			}
			minAngles[i] =  minAngle;
			minIndexs[i] = index[1];
		}
		while(size > solutionSize){
			//System.out.println("index0= "+index[0]+", index1= "+index[1]);
			double[] lamda = getMiddleAngleVector(solutionSet.get(index[0]),solutionSet.get(index[1]));
			double fitness0 = comptuePBIFitness(solutionSet.get(index[0]),lamda);
			double fitness1 = comptuePBIFitness(solutionSet.get(index[1]),lamda);
			if(fitness0 > fitness1){
			    solutionSet.get(index[0]).setRemove(true);
			    /*
			     * 更新把index[0]个体当作角度最近个体i的最近个体indexs[i];
			     */
			    for(int a=0;a<solutionSet.size();a++){
			    	if((minIndexs[a] == index[0]) && !solutionSet.get(a).isRemove()){
			    		double min = Double.MAX_VALUE;
			    		int sb = -1;
			    		double ss = 0.0;
				    	for(int b=0;b<solutionSet.size();b++){
				    		if((!solutionSet.get(b).isRemove()) && (b!=a)){
				    		ss = computeAngle(solutionSet.get(a),solutionSet.get(b));
				    			if(min > ss){
				    				min = ss;
				    				sb = b;
				    			}
				    		}
				    	}
				    	minAngles[a] = min;
				    	minIndexs[a] = sb;
				    }
			    }
			}else{
				solutionSet.get(index[1]).setRemove(true);
				/*
			     * 更新把index[1]个体当作角度最近个体i的最近个体indexs[i];
			     */
				   for(int c=0;c<solutionSet.size();c++){
				    	if((minIndexs[c] == index[1]) && !solutionSet.get(c).isRemove()){
				    		double min = Double.MAX_VALUE;
				    		int sb = -1;
				    		double ss = 0.0;
					    	for(int d=0;d<solutionSet.size();d++){
					    		if((!solutionSet.get(d).isRemove()) && (d != c)){
					    			ss = computeAngle(solutionSet.get(c),solutionSet.get(d));
					    			if(min > ss){
					    				min = ss;
					    				sb = d;
					    			}
					    		}
					    	}
					    	minAngles[c] = min;
					    	minIndexs[c] = sb;
					    }
				    }
			}
			/*
		     * 更新当前最近角度的两个个体的index;
		     */
			double sAngle = Double.MAX_VALUE;
			for(int p=0;p<solutionSet.size();p++){
				if(!solutionSet.get(p).isRemove()){
					if(sAngle > minAngles[p]){
						sAngle = minAngles[p];
						index[0] = p;
						index[1] = minIndexs[p]; 	
					}
				}
			}
			size--;
		}
	}
    
	public double[] getMiddleAngleVector(Solution s1,Solution s2){
		int objectiveSize = s1.getNumberOfObjectives();
		double[] lamda = new double[objectiveSize];
		double d1 = s1.getDistanceToIdealPoint();//S1到理想点的距离
		double d2 = s2.getDistanceToIdealPoint();//S2到理想点的距离
		for(int i=0; i<objectiveSize; i++){
			lamda[i] = s1.getNormalizedObjective(i)/d1 + s2.getNormalizedObjective(i)/d2;
		}
		return lamda;
	}
	
	public double comptuePBIFitness(Solution sol,double[] lamda){
		double fitness = 0.0;
		double d1, d2, norm;

		d1 = d2 = norm = 0.0;
		for (int i = 0; i < sol.numberOfObjectives(); i++) {
			d1 += (sol.getNormalizedObjective(i)) * lamda[i];
			norm += lamda[i]*lamda[i];
		}
		norm = Math.sqrt(norm);
		d1 = Math.abs(d1)/norm;
		for(int j=0; j<sol.numberOfObjectives(); j++){
			d2 += (sol.getNormalizedObjective(j)-d1*(lamda[j]/norm))
					*(sol.getNormalizedObjective(j)-d1*(lamda[j]/norm));
		}
		d2 = Math.sqrt(d2);
		fitness = d1 + 10*d2;
		return fitness;
	}
	
	 public void addSolutionBasedAngle(SolutionSet population,SolutionSet lastFront,int remainSize){
	    	int remain = remainSize;
	    	int obj = this.problem_.getNumberOfObjectives();
	    	int n = lastFront.size();   
			double [] angles = new double[n];
			int [] index = new int[n];
			boolean [] removed = new boolean[n];
			for(int r=0;r<n;r++){
				removed[r] = false;
			}
			
			/**
			  *  If the population is empty that is common in higher objective space
			  */
			if (population.size()==0) {
				
				for (int o = 0; o < obj;o++) {
					double minAngle2Axis = 1.0e+30;
					int minAngle2AxisID = -1;
					
					for (int i = 0;i < n;i++) {										
						if (removed[i] == false ) {
							Solution solLastFront = lastFront.get(i);
							double angle = Math.acos(Math.abs(solLastFront.getNormalizedObjective(o)/solLastFront.getDistanceToIdealPoint()) );
						
							if (angle < minAngle2Axis) {
								minAngle2Axis = angle;
								minAngle2AxisID = i;
							}
						}	
					}// for 
					 population.add(new Solution(lastFront.get(minAngle2AxisID)));
					 removed[minAngle2AxisID] = true;
					 remain--;
				} // for o 						
			} // If population.size == 0
			/**
			 * Associate each solution in the last front with a solution in the population
			 */
			for (int i = 0;i < n;i++) {
				if(!removed[i]){
					Solution solLastFront = lastFront.get(i);
					double minAng =  1.0e+30;
					int minAngID = -1;
					
					for (int j=0;j < population.size();j++){
						Solution solPop = population.get(j);
						double angle = computeAngle(solLastFront,solPop);
						if (angle < minAng){
							minAng = angle;
							minAngID = j;						
						}					
					} // for j		
					angles[i] = minAng;	
					index[i] = minAngID;
				}
			} // for i	
			
			/**
			 * Niching procedure
			 */
			while(remain > 0){
				/**
				 * Step 1: Find max and min angles 
				 */
				 int maxAngleID = -1;
				 double maxAngle = -1.0e+30;
				 
				 int minAglID = -1;
				 double minAgl = 1.0e+30;
				 
				 for(int j = 0;j < n;j++){
					// Find max angle 
					 if (!removed[j] && angles[j] > maxAngle){
						 maxAngle = angles[j];
						 maxAngleID = j;
					 }
					 
					// Find min angle 
					 if (!removed[j] && angles[j] < minAgl){
						 minAgl = angles[j];
						 minAglID = j;
					 }
				 } // for
				 
				 /**
				  * Step 2: Maximum-angle-first principle 			  
				  * 				  
				 */
					 
				 if (maxAngleID != -1 && !removed[maxAngleID]) {// Not all solutions in the last front have been added
					 population.add(lastFront.get(maxAngleID));
					 removed[maxAngleID] = true;
					 remain--;
					 
					 // Update angles			 
					 for (int i = 0;i < n; i++){ // For each solution in the last front
						 if (!removed[i]) {
							 double angle = computeAngle(lastFront.get(i),lastFront.get(maxAngleID));
							 if (angle < angles[i]){
								 angles[i] = angle;
								 index[i] = population.size() - 1;
							 }
						 }//if				 
					 }//for i	
					 
				 } else {
					 System.out.println("All solutions in the last front have been added!");
					 break;
				 }	 
			}
			
		}
	 
	 void getNextPopulation(SolutionSet pop){
			Ranking ranking = new ThetaRanking(pop, lambda_, zideal_, 
					2.0, true);
			
			int remain = populationSize_;
			int index = 0;
			SolutionSet front = null;
			population_.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);

			
			while ((remain > 0) && (remain >= front.size())) {
				
				
				for (int k = 0; k < front.size(); k++) {
					population_.add(front.get(k));
				} // for

				// Decrement remain
				remain = remain - front.size();

				// Obtain the next front
				index++;
				if (remain > 0) {
					front = ranking.getSubfront(index);
				} // if
			} // while

			if (remain > 0) { // front contains individuals to insert
			
				int[] perm = new Permutation().intPermutation(front.size());
				for (int k = 0; k < remain; k++) {
					population_.add(front.get(perm[k]));
				} // for
				remain = 0;
				
			} // if
		}
	 
	  SolutionSet[] getParetoFronts() {
			
			SolutionSet[] sets = new SolutionSet[2];
			Ranking ranking = new NondominatedRanking(union_);

			int remain = populationSize_;
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

}
