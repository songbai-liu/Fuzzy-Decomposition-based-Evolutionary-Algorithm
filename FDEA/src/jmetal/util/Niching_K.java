package jmetal.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import Jama.Matrix;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;

public class Niching_K {
	SolutionSet population;
	SolutionSet lastFront;
	SolutionSet mgPopulation;

	SolutionSet union;

	int obj;
	int remain;
	
	double[][] lambda;

	public Niching_K(SolutionSet population, SolutionSet lastFront,
			double[][] lambda, int remain) {

		this.population = population;
		this.lastFront = lastFront;

		this.remain = remain;
		this.lambda = lambda;

		this.mgPopulation = population.union(lastFront);

		if (population.size() > 0)
			this.obj = population.get(0).numberOfObjectives();
		else
			this.obj = lastFront.get(0).numberOfObjectives();
	}

	public void execute() {	
		associate();
		assignment();
	}

	

	public void associate() {

		for (int k = 0; k < mgPopulation.size(); k++) {

			Solution sol = mgPopulation.get(k);

			double min = calVDistance(sol, lambda[0]);
			int index = 0;

			for (int j = 1; j < lambda.length; j++) {
				double dist = calVDistance(sol, lambda[j]);
				if (dist < min) {
					min = dist;
					index = j;
				}
			}
			sol.setClusterID(index);
			sol.setVDistance(min);
		}

	}

	public void assignment() {
		int[] ro = new int[lambda.length];
		boolean[] flag = new boolean[lambda.length];

		for (int k = 0; k < population.size(); k++) {
			ro[population.get(k).getClusterID()]++;
		}

		int num = 0;

		while (num < remain) {
			int[] perm = new Permutation().intPermutation(ro.length);

			int min = Integer.MAX_VALUE;
			int id = -1;

			for (int i = 0; i < perm.length; i++) {
				if ((!flag[perm[i]]) && (ro[perm[i]] < min)) {
					min = ro[perm[i]];
					id = perm[i];
				}
			}

			List<Integer> list = new ArrayList<Integer>();

			for (int k = 0; k < lastFront.size(); k++) {
				if (lastFront.get(k).getClusterID() == id)
					list.add(k);
			}

			if (list.size() != 0) {
				int index = 0;
				if (ro[id] == 0) {
					double minDist = Double.MAX_VALUE;

					for (int j = 0; j < list.size(); j++) {
						if (lastFront.get(list.get(j)).getVDistance() < minDist) {
							minDist = lastFront.get(list.get(j)).getVDistance();
							index = j;
						}
					}
				} else {
					index = PseudoRandom.randInt(0, list.size() - 1);
				}

				population.add(lastFront.get(list.get(index)));
				ro[id]++;

				lastFront.remove(list.get(index));
				num++;
			} else {
				flag[id] = true;
			}

		}
	}

	
	public double calVDistance(Solution sol, double[] ref){
		return calNormlizedVDistance(sol, ref);
	}
	
	public double calNormlizedVDistance(Solution sol, double[] ref) {

		double ip = 0;
		double refLenSQ = 0;

		for (int j = 0; j < obj; j++) {

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
}
