//  Ranking.java
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

package jmetal.util;

import jmetal.core.SolutionSet;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.comparators.OverallConstraintViolationComparator;
import jmetal.util.comparators.ReducedDominanceComparator;

import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * This class implements some facilities for ranking solutions. Given a
 * <code>SolutionSet</code> object, their solutions are ranked according to
 * scheme proposed in NSGA-II; as a result, a set of subsets are obtained. The
 * subsets are numbered starting from 0 (in NSGA-II, the numbering starts from
 * 1); thus, subset 0 contains the non-dominated solutions, subset 1 contains
 * the non-dominated solutions after removing those belonging to subset 0, and
 * so on.
 */
public class ReducedRanking {

	/**
	 * The <code>SolutionSet</code> to rank
	 */
	private SolutionSet solutionSet_;

	/**
	 * An array containing all the fronts found during the search
	 */
	private SolutionSet[] ranking_;

	/**
	 * stores a <code>Comparator</code> for dominance checking
	 */
	private static final Comparator dominance_ = new ReducedDominanceComparator();

	/**
	 * stores a <code>Comparator</code> for Overal Constraint Violation
	 * Comparator checking
	 */
	private static final Comparator constraint_ = new OverallConstraintViolationComparator();
	
	//private DominanceComparator1 dominance1_ = new DominanceComparator1(group_);

	/**
	 * Constructor.
	 * 
	 * @param solutionSet
	 *            The <code>SolutionSet</code> to be ranked.
	 */
	public ReducedRanking(SolutionSet solutionSet) {
		for(int p=0;p<solutionSet.size();p++){
    		solutionSet.get(p).setRemove(false);
    	}
		//solutionSet_ = new SolutionSet(solutionSet.size());
		solutionSet_ = solutionSet;
	/*	for(int p=0;p<solutionSet.size();p++){
			solutionSet_.add(solutionSet.get(p));
		}*/
		sortSolutionBsedSumValue(solutionSet_);
		
		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[solutionSet_.size() + 1];
		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();
		
		// flagDominate is an auxiliar encodings.variable
		int id = 0;
		//ranking_[id].add(solutionSet_.get(0));
		//solutionSet_.get(0).setRank(0);
		//front[id].add(0);
		//solutionSet_.remove(0);
		//solutionSet_.get(0).setRemove(true);
		int size = solutionSet_.size();
		while(size > 0){
			for(int p=0;p<solutionSet_.size();p++){
				if(!solutionSet_.get(p).isRemove()){
					front[id].add(p);
					solutionSet_.get(p).setRank(id);
					solutionSet_.get(p).setRemove(true);
					size--;
					break;
				}
			}
			for(int i=0;i<solutionSet_.size();i++){
				
				//for(int j=0;j<ranking_[id].size();j++){
				if(!solutionSet_.get(i).isRemove()){
					int flagDominate = -2;
					for(int j=0;j<front[id].size();j++){
						flagDominate = dominance_.compare(solutionSet_.get(i),
								solutionSet_.get(front[id].get(j)));
						if(flagDominate == 1){
							break;
						}
					}
					if(flagDominate != 1){
						//ranking_[id].add(solutionSet_.get(i));
						front[id].add(i);
						solutionSet_.get(i).setRank(id);
						solutionSet_.get(i).setRemove(true);
						size--;
					}
				}
			}
			id++;
		}
		Iterator<Integer> iterator; // Iterators
		ranking_ = new SolutionSet[id];
		for(int s=0;s<id;s++){
			ranking_[s] = new SolutionSet(front[s].size());
			iterator = front[s].iterator();
			while (iterator.hasNext()) {
				ranking_[s].add(solutionSet_.get(iterator.next()));
			}
		}
		
	} // Ranking
	
	public void sortSolutionBsedSumValue(SolutionSet solutionSet){
		for(int i=0;i<solutionSet.size();i++){
			double sumValue = 0.0;
			for(int j=0;j<solutionSet.get(0).numberOfObjectives();j++){
				sumValue += solutionSet.get(i).getNormalizedObjective(j);
			}
			solutionSet.get(i).setFitness(sumValue);
		}
		solutionSet.sort(new FitnessComparator());
	}
	
	/**
	 * Returns a <code>SolutionSet</code> containing the solutions of a given
	 * rank.
	 * 
	 * @param rank
	 *            The rank
	 * @return Object representing the <code>SolutionSet</code>.
	 */
	public SolutionSet getSubfront(int rank) {
		return ranking_[rank];
	} // getSubFront

	/**
	 * Returns the total number of subFronts founds.
	 */
	public int getNumberOfSubfronts() {
		return ranking_.length;
	} // getNumberOfSubfronts
	
/*	public static void main(String[] args){
		List<Integer>[] front = new List[2];
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();
		front[0].add(1);
		front[0].add(1);
		front[0].add(1);
		front[1].add(1);
		front[1].add(1);
		System.out.println(front[0].size());
		System.out.println(front[1].size());
	}*/
} // Ranking
