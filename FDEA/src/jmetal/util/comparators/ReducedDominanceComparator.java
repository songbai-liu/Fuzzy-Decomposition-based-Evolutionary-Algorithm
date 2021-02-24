//  DominanceComparator.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
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

package jmetal.util.comparators;

import jmetal.core.Solution;
import jmetal.util.PseudoRandom;

import java.util.Comparator;

/**
 * This class implements a <code>Comparator</code> (a method for comparing
 * <code>Solution</code> objects) based on a constraint violation test + 
 * dominance checking, as in NSGA-II.
 */
public class ReducedDominanceComparator implements Comparator{
	//implements Comparator
  IConstraintViolationComparator violationConstraintComparator_ ;

  /**
   * Constructor
   */
  public ReducedDominanceComparator() {
    violationConstraintComparator_ = new OverallConstraintViolationComparator(); 
    //violationConstraintComparator_ = new NumberOfViolatedConstraintComparator(); 
  }

  /**
   * Constructor
   * @param comparator
   */
  public ReducedDominanceComparator(IConstraintViolationComparator comparator) {
    violationConstraintComparator_ = comparator ;
  }
 
 /**
  * Compares two solutions.
  * @param object1 Object representing the first <code>Solution</code>.
  * @param object2 Object representing the second <code>Solution</code>.
  * @return -1, or 0, or 1 if solution1 dominates solution2, both are 
  * non-dominated, or solution1  is dominated by solution22, respectively.
  */
  public int compare(Object object1, Object object2) {
    if (object1==null)
      return 1;
    else if (object2 == null)
      return -1;
    
    Solution solution1 = (Solution)object1;
    Solution solution2 = (Solution)object2;

    int dominate1 ; // dominate1 indicates if some objective of solution1 
                    // dominates the same objective in solution2. dominate2
    int dominate2 ; // is the complementary of dominate1.

    dominate1 = 0 ; 
    dominate2 = 0 ;
    
    int flag; //stores the result of the comparison

    // Test to determine whether at least a solution violates some constraint
    if (violationConstraintComparator_.needToCompare(solution1, solution2))
      return violationConstraintComparator_.compare(solution1, solution2) ;
    
    if(solution1.numberOfObjectives() < 4){
    	System.out.println("It is not a many objective optimization problem!");
    	System.exit(0);
    }
    //random selection 3 objective index for compare
    int[] id = new int[3];
    id[0] = (int) (PseudoRandom.randInt(0,solution1.numberOfObjectives() - 1));
    do{
    	id[1] = (int) (PseudoRandom.randInt(0,solution1.numberOfObjectives() - 1));
    }while(id[0] == id[1]);
    do{
    	id[2] = (int) (PseudoRandom.randInt(0,solution1.numberOfObjectives() - 1));
    }while(id[0] == id[2] || id[2] == id[1]);
    
    // Equal number of violated constraints. Applying a dominance Test then
    double value1, value2;
    for (int i = 0; i < 3; i++) {
    	double rd = PseudoRandom.randDouble(0.3, 0.7);
    	value1 = rd*solution1.getNormalizedObjective(id[i]);
        value2 = rd*solution2.getNormalizedObjective(id[i]);
    	for(int j=0;j<solution1.numberOfObjectives();j++){
    		if(j != id[i]){
    			value1 += ((1-rd)/(solution1.numberOfObjectives()-1))*solution1.getNormalizedObjective(j);
    			value2 += ((1-rd)/(solution1.numberOfObjectives()-1))*solution2.getNormalizedObjective(j);
    		}
    	}
      //value1 = solution1.getObjective(id[i]);
      //value2 = solution2.getObjective(id[i]);
      if (value1 < value2) {
        flag = -1;
      } else if (value1 > value2) {
        flag = 1;
      } else {
        flag = 0;
      }
      
      if (flag == -1) {
        dominate1 = 1;
      }
      
      if (flag == 1) {
        dominate2 = 1;           
      }
    }
            
    if (dominate1 == dominate2) {            
      return 0; //No one dominate the other
    }
    if (dominate1 == 1) {
      return -1; // solution1 dominate
    }
    return 1;    // solution2 dominate   
  } // compare  
} // DominanceComparator
