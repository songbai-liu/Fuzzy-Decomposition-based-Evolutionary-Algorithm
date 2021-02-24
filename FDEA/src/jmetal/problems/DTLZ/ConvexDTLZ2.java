//  DTLZ2Convex.java
//
//  Author:
//       Yi Xiang <gzhuxiang_yi@163.com>



package jmetal.problems.DTLZ;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem ZDT2
 */
public class ConvexDTLZ2 extends Problem{
       
 /** 
  * Creates a default convex DTLZ2  problem (12 variables and 3 objectives)
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public ConvexDTLZ2(String solutionType) throws ClassNotFoundException {
    this(solutionType, 12, 3);//3 objectives 
//    this(solutionType, 13, 4);//4 objectives 
//    this(solutionType, 14, 5);//5 objectives 
//    this(solutionType, 15, 6);//6 objectives 
//    this(solutionType, 17, 8);//8 objectives 
//    this(solutionType, 19, 10);//10 objectives 
//    this(solutionType, 24, 15);//15 objectives  
  
  } // DTLZ2Convex

 /**
  * Creates a new instance of convex DTLZ2 
  * @param numberOfVariables Number of variables
  * @param numberOfObjectives Number of objective functions
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public ConvexDTLZ2(String  solutionType,
               Integer numberOfVariables,
               Integer numberOfObjectives) {
    numberOfVariables_  = numberOfVariables;
    numberOfObjectives_ = numberOfObjectives;
    numberOfConstraints_= 0;
    problemName_        = "ConvexDTLZ2";
        
    upperLimit_ = new double[numberOfVariables_];
    lowerLimit_ = new double[numberOfVariables_];        
    for (int var = 0; var < numberOfVariables_; var++){
      lowerLimit_[var] = 0.0;
      upperLimit_[var] = 1.0;
    } //for
        
    if (solutionType.compareTo("BinaryReal") == 0)
    	solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compareTo("Real") == 0)
    	solutionType_ = new RealSolutionType(this) ;
    else {
    	System.out.println("Error: solution type " + solutionType + " invalid") ;
    	System.exit(-1) ;
    }            
  } //DTLZ2Convex
        
  /** 
  * Evaluates a solution 
  * @param solution The solution to evaluate
   * @throws JMException 
  */    
  public void evaluate(Solution solution) throws JMException {
    Variable[] gen  = solution.getDecisionVariables();

    double [] x = new double[numberOfVariables_];
    double [] f = new double[numberOfObjectives_];
    int k = numberOfVariables_ - numberOfObjectives_ + 1;
        
    for (int i = 0; i < numberOfVariables_; i++)
      x[i] = gen[i].getValue();
        
    double g = 0.0;
    for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
      g += (x[i] - 0.5)*(x[i] - 0.5);
        
    for (int i = 0; i < numberOfObjectives_; i++)
      f[i] = 1.0 + g;
        
    for (int i = 0; i < numberOfObjectives_; i++){
      for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)            
        f[i] *= Math.cos(x[j]*0.5*Math.PI);                
        if (i != 0){
          int aux = numberOfObjectives_ - (i + 1);
          f[i] *= Math.sin(x[aux]*0.5*Math.PI);
        } //if 
    } // for
        
    // Map objectives as fi = fi^4,i=1,..,M-1, fM = fM^2
    for (int i = 0; i < numberOfObjectives_ - 1; i++) {
    	f[i] = Math.pow(f[i], 4);
    }
    
    f[numberOfObjectives_-1] = Math.pow(f[numberOfObjectives_-1], 2);	
    	
    for (int i = 0; i < numberOfObjectives_; i++)
      solution.setObjective(i,f[i]);        
  }    
} //evaluate
