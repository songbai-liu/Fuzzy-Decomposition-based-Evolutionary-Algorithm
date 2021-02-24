//  Solution.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Description: 
// 
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.core;

import jmetal.encodings.variable.Binary;

import java.io.Serializable;

/**
 * Class representing a solution for a problem.
 */
public class Solution implements Serializable {
	/**
	 * Stores the problem
	 */
	private Problem problem_;

	/**
	 * Stores the type of the encodings.variable
	 */
	private SolutionType type_;

	/**
	 * Stores the decision variables of the solution.
	 */
	private Variable[] variable_;

	/**
	 * Stores the objectives values of the solution.
	 */
	private final double[] objective_;
	
	private double[] translatedObjectives;

	/**
	 * Stores the number of objective values of the solution
	 */
	private int numberOfObjectives_;

	public int getNumberOfObjectives() {
		return numberOfObjectives_;
	}
	/**
	 * Stores the so called fitness value. Used in some metaheuristics
	 */
	private double fitness_;

	/**
	 * Used in algorithm AbYSS, this field is intended to be used to know when a
	 * <code>Solution</code> is marked.
	 */
	private boolean marked_;

	/**
	 * Stores the so called rank of the solution. Used in NSGA-II
	 */
	private int rank_;

	/**
	 * Stores the overall constraint violation of the solution.
	 */
	private double overallConstraintViolation_;

	/**
	 * Stores the number of constraints violated by the solution.
	 */
	private int numberOfViolatedConstraints_;

	/**
	 * This field is intended to be used to know the location of a solution into
	 * a <code>SolutionSet</code>. Used in MOCell
	 */
	private int location_;
	
	private double diversity_;
	
	private double associateDist_;

	/**
	 * This field is intended to be used to know the region of a solution
	 * <code>SolutionSet</code>. Used in MST
	 */
	private int region_;

	/**
	 * Stores the distance to his k-nearest neighbor into a
	 * <code>SolutionSet</code>. Used in SPEA2.
	 */
	private double kDistance_;

	/**
	 * Stores the crowding distance of the the solution in a
	 * <code>SolutionSet</code>. Used in NSGA-II.
	 */
	private double crowdingDistance_;
	
	private int localDensity_;
	private double gaussianLocalDensity_;
	private double centerDistance_;

	/**
	 * Stores the distance between this solution and a <code>SolutionSet</code>.
	 * Used in AbySS.
	 */
	private double distanceToSolutionSet_;

	private double max_distance_;
	  
	private double min_distance_;
	
	private int cross_type;
	
	private int clone_num;
	
	private int clusterID_;
	private double vDistance_;
	
	private Solution lamda_;
	
	public Solution getLamda_() {
		return lamda_;
	}
	public void setLamda_(Solution lamda_) {
		this.lamda_ = lamda_;
	}
	private double[] normalizedObjective_;
	private double[] unitHyperplaneObjective_;
	
	private double[] unitHypersphereObjective_;
	
	public double getUnitHyperplaneObjective(int i) {
		return unitHyperplaneObjective_[i];
	}
	public void setUnitHyperplaneObjective(int i,double unitValue) {
		unitHyperplaneObjective_[i] = unitValue;
	}
	
	public double getUnitHypersphereObjective(int i) {
		return unitHypersphereObjective_[i];
	}
	public void setUnitHypersphereObjective(int i,double unitValue) {
		unitHypersphereObjective_[i] = unitValue;
	}
	/*.........Used in NMPSO.............*/
	private double convergenceDistance_;
	public double getConvergenceDistance_() {
		return convergenceDistance_;
	}
	public void setConvergenceDistance_(double convergenceDistance_) {
		this.convergenceDistance_ = convergenceDistance_;
	}
	/*.....................distance to ideal point......................*/	
	private double distanceToIdealPoint;
	public double getDistanceToIdealPoint() {
		return distanceToIdealPoint;
	}
	public void setDistanceToIdealPoint(double distanceToIdealPoint) {
		this.distanceToIdealPoint = distanceToIdealPoint;
	}
	/*..........................end.................................*/
	
	/*.....................Sum Value......................*/	
	private double sumValue;
	public double getSumValue() {
		return sumValue;
	}
	public void setSumValue(double sumValue) {
		this.sumValue = sumValue;
	}
	/*..........................end.................................*/
	
	/*.....................distance to nadir point......................*/
	
	private double distanceToNadirPoint;
	public double getDistanceToNadirPoint() {
		return distanceToNadirPoint;
	}
	public void setDistanceToNadirPoint(double distanceToNadirPoint) {
		this.distanceToNadirPoint = distanceToNadirPoint;
	}
	/*..........................end................................*/
	
	/*...................if the solution removed...................*/
	private boolean remove;
	public boolean isRemove() {
		return remove;
	}
	public void setRemove(boolean remove) {
		this.remove = remove;
	}
	/*.............................................................*/
	
	/*...................if the solution removed...................*/
	private boolean inDiversitySet;
	private boolean inConvergenceSet;
	public boolean isInDiversitySet() {
		return this.inDiversitySet;
	}
	public void setInDiversitySet(boolean id1) {
		this.inDiversitySet = id1;
	}
	
	public boolean isInConvergenceSet() {
		return this.inConvergenceSet;
	}
	public void setInConvergenceSet(boolean id2) {
		this.inConvergenceSet = id2;
	}
	/*.............................................................*/
	/**
	 * Constructor.
	 */
	public Solution() {
		problem_ = null;
		marked_ = false;
		remove = false;
		overallConstraintViolation_ = 0.0;
		numberOfViolatedConstraints_ = 0;
		type_ = null;
		variable_ = null;
		objective_ = null;
		this.translatedObjectives = null;
		lamda_ = this;
		inDiversitySet = false;
		inConvergenceSet = false;
	} // Solution

	/**
	 * Constructor
	 * 
	 * @param numberOfObjectives
	 *            Number of objectives of the solution
	 * 
	 *            This constructor is used mainly to read objective values from
	 *            a file to variables of a SolutionSet to apply quality
	 *            indicators
	 */
	public Solution(int numberOfObjectives) {
		numberOfObjectives_ = numberOfObjectives;
		objective_ = new double[numberOfObjectives];
		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		unitHypersphereObjective_ = new double[numberOfObjectives_];
		this.translatedObjectives = new double[this.numberOfObjectives_];
		lamda_ = this;
		inDiversitySet = false;
		inConvergenceSet = false;
	}

	/**
	 * Constructor.
	 * 
	 * @param problem
	 *            The problem to solve
	 * @throws ClassNotFoundException
	 */
	public Solution(Problem problem) throws ClassNotFoundException {
		problem_ = problem;
		type_ = problem.getSolutionType();
		numberOfObjectives_ = problem.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		//clone_num=new int[numberOfObjectives_];
		
		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		unitHypersphereObjective_ = new double[numberOfObjectives_];
		remove = false;
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;
		// <-
		// variable_ = problem.solutionType_.createVariables() ;
		variable_ = type_.createVariables();
		this.translatedObjectives = new double[this.numberOfObjectives_];
		lamda_ = this;
		inDiversitySet = false;
		inConvergenceSet = false;
	} // Solution
	
	public Solution(Problem problem,int groupSize) throws ClassNotFoundException {
		problem_ = problem;
		type_ = problem.getSolutionType();
		numberOfObjectives_ = problem.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		//clone_num=new int[numberOfObjectives_];
		
		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		unitHypersphereObjective_ = new double[numberOfObjectives_];
		remove = false;
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;
		// variable_ = problem.solutionType_.createVariables() ;
		variable_ = type_.createVariables();
		lamda_ = this;
		inDiversitySet = false;
		inConvergenceSet = false;
		
	} // Solution

	static public Solution getNewSolution(Problem problem)
			throws ClassNotFoundException {
		return new Solution(problem);
	}

	/**
	 * Constructor
	 * 
	 * @param problem
	 *            The problem to solve
	 */
	public Solution(Problem problem, Variable[] variables) {
		problem_ = problem;
		type_ = problem.getSolutionType();
		numberOfObjectives_ = problem.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		
		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		unitHypersphereObjective_ = new double[numberOfObjectives_];
		
		remove = false;
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;

		variable_ = variables;
		this.translatedObjectives = new double[this.numberOfObjectives_];
		lamda_ = this;
		inDiversitySet = false;
		inConvergenceSet = false;
	} // Constructor

	/**
	 * Copy constructor.
	 * 
	 * @param solution
	 *            Solution to copy.
	 */
	public Solution(Solution solution) {
		problem_ = solution.problem_;
		type_ = solution.type_;

		numberOfObjectives_ = solution.numberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		
		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		unitHypersphereObjective_ = new double[numberOfObjectives_];
		
		for (int i = 0; i < objective_.length; i++) {
			objective_[i] = solution.getObjective(i);
			normalizedObjective_[i] = solution.getNormalizedObjective(i);
			unitHyperplaneObjective_[i] = solution.getUnitHyperplaneObjective(i);
			unitHypersphereObjective_[i] = solution.getUnitHypersphereObjective(i);
			//this.translatedObjectives[i] = solution.getIthTranslatedObjective(i);
		} // for
			// <-
		this.translatedObjectives = new double[this.numberOfObjectives_];
		remove = false;
		variable_ = type_.copyVariables(solution.variable_);
		overallConstraintViolation_ = solution.getOverallConstraintViolation();
		numberOfViolatedConstraints_ = solution.getNumberOfViolatedConstraint();
		distanceToSolutionSet_ = solution.getDistanceToSolutionSet();
		distanceToIdealPoint = solution.getDistanceToIdealPoint();
		distanceToNadirPoint = solution.getDistanceToNadirPoint();
		sumValue = solution.getSumValue();
		crowdingDistance_ = solution.getCrowdingDistance();
		localDensity_ = solution.getLocalDensity();
		gaussianLocalDensity_ = solution.getGaussianLocalDensity();
		centerDistance_ = solution.getCenterDistance();
		convergenceDistance_ = solution.getConvergenceDistance_();
		kDistance_ = solution.getKDistance();
		fitness_ = solution.getFitness();
		marked_ = solution.isMarked();
		rank_ = solution.getRank();
		location_ = solution.getLocation();
		max_distance_=  solution.getmaxDistance();
		min_distance_=  solution.getminDistance();
		cross_type=  solution.getcross_type();
		clone_num=  solution.getclone_num();
		lamda_ = this;
		
		inDiversitySet = false;
		inConvergenceSet = false;
	} // Solution
	
	public double getGaussianLocalDensity() {
		return gaussianLocalDensity_;
	}
	public void setGaussianLocalDensity(double gaussianLocalDensity_) {
		this.gaussianLocalDensity_ = gaussianLocalDensity_;
	}
	public double getCenterDistance() {
		return centerDistance_;
	}
	public void setCenterDistance(double centerDistance_) {
		this.centerDistance_ = centerDistance_;
	}
	public int getLocalDensity() {
		return localDensity_;
	}
	public void setLocalDensity(int localDensity_) {
		this.localDensity_ = localDensity_;
	}
	public double getCrowdingDistance_() {
		return crowdingDistance_;
	}

	public void setCrowdingDistance_(double crowdingDistance_) {
		this.crowdingDistance_ = crowdingDistance_;
	}

	/**
	 * Sets the distance between this solution and a <code>SolutionSet</code>.
	 * The value is stored in <code>distanceToSolutionSet_</code>.
	 * 
	 * @param distance
	 *            The distance to a solutionSet.
	 */
	public void setDistanceToSolutionSet(double distance) {
		distanceToSolutionSet_ = distance;
	} // SetDistanceToSolutionSet

	public double getmaxDistance(){
	    return max_distance_;
	  } // getDistanceToSolutionSet
	public void setmaxDistance(double distance){
		    max_distance_ = distance;
		  } // SetDistanceToSolutionSet

	public double getminDistance(){
		    return min_distance_;
		  } // getDistanceToSolutionSet
	public void setminDistance(double distance){
	    min_distance_ = distance;
    } // SetDistanceToSolutionSet
		  
	 public int getcross_type(){
				 return cross_type;
			} // getDistanceToSolutionSet
	public void setcross_type(int type){
				cross_type = type;
			} // SetDistanceToSolutionSet  
			
	public int getclone_num(){
				 return clone_num;
			} // getDistanceToSolutionSet
			public void setclone_num(int num){
				clone_num = num;
			} // SetDistanceToSolutionSet  
	/**
	 * Gets the distance from the solution to a <code>SolutionSet</code>. <b>
	 * REQUIRE </b>: this method has to be invoked after calling
	 * <code>setDistanceToPopulation</code>.
	 * 
	 * @return the distance to a specific solutionSet.
	 */
	public double getDistanceToSolutionSet() {
		return distanceToSolutionSet_;
	} // getDistanceToSolutionSet

	/**
	 * Sets the distance between the solution and its k-nearest neighbor in a
	 * <code>SolutionSet</code>. The value is stored in <code>kDistance_</code>.
	 * 
	 * @param distance
	 *            The distance to the k-nearest neighbor.
	 */
	public void setKDistance(double distance) {
		kDistance_ = distance;
	} // setKDistance

	/**
	 * Gets the distance from the solution to his k-nearest nighbor in a
	 * <code>SolutionSet</code>. Returns the value stored in
	 * <code>kDistance_</code>. <b> REQUIRE </b>: this method has to be invoked
	 * after calling <code>setKDistance</code>.
	 * 
	 * @return the distance to k-nearest neighbor.
	 */
	public double getKDistance() {
		return kDistance_;
	} // getKDistance

	/**
	 * Sets the crowding distance of a solution in a <code>SolutionSet</code>.
	 * The value is stored in <code>crowdingDistance_</code>.
	 * 
	 * @param distance
	 *            The crowding distance of the solution.
	 */
	public void setCrowdingDistance(double distance) {
		crowdingDistance_ = distance;
	} // setCrowdingDistance

	/**
	 * Gets the crowding distance of the solution into a
	 * <code>SolutionSet</code>. Returns the value stored in
	 * <code>crowdingDistance_</code>. <b> REQUIRE </b>: this method has to be
	 * invoked after calling <code>setCrowdingDistance</code>.
	 * 
	 * @return the distance crowding distance of the solution.
	 */
	public double getCrowdingDistance() {
		return crowdingDistance_;
	} // getCrowdingDistance

	/**
	 * Sets the fitness of a solution. The value is stored in
	 * <code>fitness_</code>.
	 * 
	 * @param fitness
	 *            The fitness of the solution.
	 */
	public void setFitness(double fitness) {
		fitness_ = fitness;
	} // setFitness

	/**
	 * Gets the fitness of the solution. Returns the value of stored in the
	 * encodings.variable <code>fitness_</code>. <b> REQUIRE </b>: This method
	 * has to be invoked after calling <code>setFitness()</code>.
	 * 
	 * @return the fitness.
	 */
	public double getFitness() {
		return fitness_;
	} // getFitness

	/**
	 * Sets the value of the i-th objective.
	 * 
	 * @param i
	 *            The number identifying the objective.
	 * @param value
	 *            The value to be stored.
	 */
	public void setObjective(int i, double value) {
		objective_[i] = value;
	} // setObjective

	/**
	 * Returns the value of the i-th objective.
	 * 
	 * @param i
	 *            The value of the objective.
	 */
	public double getObjective(int i) {
		return objective_[i];
	} // getObjective

	/**
	 * Returns the number of objectives.
	 * 
	 * @return The number of objectives.
	 */
	public int numberOfObjectives() {
		if (objective_ == null)
			return 0;
		else
			return numberOfObjectives_;
	} // numberOfObjectives

	/**
	 * Returns the number of decision variables of the solution.
	 * 
	 * @return The number of decision variables.
	 */
	public int numberOfVariables() {
		return problem_.getNumberOfVariables();
	} // numberOfVariables

	/**
	 * Returns a string representing the solution.
	 * 
	 * @return The string.
	 */
	public String toString() {
		String aux = "";
		for (int i = 0; i < this.numberOfObjectives_; i++)
			aux = aux + this.getObjective(i) + " ";

		return aux;
	} // toString

	/**
	 * Returns the decision variables of the solution.
	 * 
	 * @return the <code>DecisionVariables</code> object representing the
	 *         decision variables of the solution.
	 */
	public Variable[] getDecisionVariables() {
		return variable_;
	} // getDecisionVariables

	/**
	 * Sets the decision variables for the solution.
	 * 
	 * @param variables
	 *            The <code>DecisionVariables</code> object representing the
	 *            decision variables of the solution.
	 */
	public void setDecisionVariables(Variable[] variables) {
		variable_ = variables;
	} // setDecisionVariables

	/**
	 * Indicates if the solution is marked.
	 * 
	 * @return true if the method <code>marked</code> has been called and, after
	 *         that, the method <code>unmarked</code> hasn't been called. False
	 *         in other case.
	 */
	public boolean isMarked() {
		return this.marked_;
	} // isMarked

	/**
	 * Establishes the solution as marked.
	 */
	public void marked() {
		this.marked_ = true;
	} // marked

	/**
	 * Established the solution as unmarked.
	 */
	public void unMarked() {
		this.marked_ = false;
	} // unMarked

	/**
	 * Sets the rank of a solution.
	 * 
	 * @param value
	 *            The rank of the solution.
	 */
	public void setRank(int value) {
		this.rank_ = value;
	} // setRank

	/**
	 * Gets the rank of the solution. <b> REQUIRE </b>: This method has to be
	 * invoked after calling <code>setRank()</code>.
	 * 
	 * @return the rank of the solution.
	 */
	public int getRank() {
		return this.rank_;
	} // getRank

	/**
	 * Sets the overall constraints violated by the solution.
	 * 
	 * @param value
	 *            The overall constraints violated by the solution.
	 */
	public void setOverallConstraintViolation(double value) {
		this.overallConstraintViolation_ = value;
	} // setOverallConstraintViolation
	
	public void setNormalizedObjective(int i, double value) {
		normalizedObjective_[i] = value;
	}

	public double getNormalizedObjective(int i) {
		return normalizedObjective_[i];
	}
	

	/**
	 * Gets the overall constraint violated by the solution. <b> REQUIRE </b>:
	 * This method has to be invoked after calling
	 * <code>overallConstraintViolation</code>.
	 * 
	 * @return the overall constraint violation by the solution.
	 */
	public double getOverallConstraintViolation() {
		return this.overallConstraintViolation_;
	} // getOverallConstraintViolation

	/**
	 * Sets the number of constraints violated by the solution.
	 * 
	 * @param value
	 *            The number of constraints violated by the solution.
	 */
	public void setNumberOfViolatedConstraint(int value) {
		this.numberOfViolatedConstraints_ = value;
	} // setNumberOfViolatedConstraint

	/**
	 * Gets the number of constraint violated by the solution. <b> REQUIRE </b>:
	 * This method has to be invoked after calling
	 * <code>setNumberOfViolatedConstraint</code>.
	 * 
	 * @return the number of constraints violated by the solution.
	 */
	public int getNumberOfViolatedConstraint() {
		return this.numberOfViolatedConstraints_;
	} // getNumberOfViolatedConstraint

	/**
	 * Sets the location of the solution into a solutionSet.
	 * 
	 * @param location
	 *            The location of the solution.
	 */
	public void setLocation(int location) {
		this.location_ = location;
	} // setLocation

	/**
	 * Gets the location of this solution in a <code>SolutionSet</code>. <b>
	 * REQUIRE </b>: This method has to be invoked after calling
	 * <code>setLocation</code>.
	 * 
	 * @return the location of the solution into a solutionSet
	 */
	public int getLocation() {
		return this.location_;
	} // getLocation
	
	
	public void setClusterID(int id){
		this.clusterID_ = id;
	}
	
	public int getClusterID(){
		return this.clusterID_;
	}
	
	public void setVDistance(double val){
		this.vDistance_ = val;
	}
	
	public double getVDistance(){
		return this.vDistance_;
	}

	/**
	 * Sets the type of the encodings.variable.
	 * 
	 * @param type
	 *            The type of the encodings.variable.
	 */
	public void setType(SolutionType type) {
		type_ = type;
	} // setType

	/**
	 * Gets the type of the encodings.variable
	 * 
	 * @return the type of the encodings.variable
	 */
	public SolutionType getType() {
		return type_;
	} // getType

	/**
	 * Returns the aggregative value of the solution
	 * 
	 * @return The aggregative value.
	 */
	public double getAggregativeValue() {
		double value = 0.0;
		for (int i = 0; i < numberOfObjectives(); i++) {
			value += getObjective(i);
		}
		return value;
	} // getAggregativeValue

	/**
	 * Returns the number of bits of the chromosome in case of using a binary
	 * representation
	 * 
	 * @return The number of bits if the case of binary variables, 0 otherwise
	 *         This method had a bug which was fixed by Rafael Olaechea
	 */
	public int getNumberOfBits() {
		int bits = 0;

		for (int i = 0; i < variable_.length; i++)
			if ((variable_[i].getVariableType() == jmetal.encodings.variable.Binary.class)
					|| (variable_[i].getVariableType() == jmetal.encodings.variable.BinaryReal.class))

				bits += ((Binary) (variable_[i])).getNumberOfBits();

		return bits;
	} // getNumberOfBits
	//-----------------Ìí¼ÓÄÚÈÝstart-MOEADD--------------------//
	public void Set_location(int i) {
		this.location_ = i;
	}

	public int read_location() {
		return this.location_;
	}

	public void Set_diversity(double i) {
		this.diversity_ = i;
	}

	public double read_diversity() {
		return this.diversity_;
	}
	
	public void setRegion(int i) {
		this.region_ = i;
	}

	public int readRegion() {
		return this.region_;
	}

	public void Set_associateDist(double distance) {
		this.associateDist_ = distance;
	}

	public double read_associateDist() {
		return this.associateDist_;
	}
		
	public void setIthTranslatedObjective(int i,double val){
		this.translatedObjectives[i] = val;
	}
	
	public double getIthTranslatedObjective(int i){
		return this.translatedObjectives[i];
	}
	
	
} // Solution
