package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;

public class TestWfgHv {
	public static void main(String[] args) throws IOException{
		double hv = 0.0;
		  int fun = 7;
		  SolutionSet population = new SolutionSet(7);
		  for(int i=0;i<7;i++){
			 Solution sol = new Solution(3);
			 population.add(sol);
		  }
		 population.get(0).setObjective(0, 0.0);
		 population.get(0).setObjective(1, 0.5);
		 population.get(0).setObjective(2, 0.5);
		 population.get(1).setObjective(0, 0.5);
		 population.get(1).setObjective(1, 0.0);
		 population.get(1).setObjective(2, 0.5);
		 population.get(2).setObjective(0, 0.5);
		 population.get(2).setObjective(1, 0.5);
		 population.get(2).setObjective(2, 0.0);
		 population.get(3).setObjective(0, 0.4);
		 population.get(3).setObjective(1, 0.4); 
		 population.get(3).setObjective(2, 0.4);
		 population.get(4).setObjective(0, 0.0);
		 population.get(4).setObjective(1, 0.0);
		 population.get(4).setObjective(2, 1.0);
		 population.get(5).setObjective(0, 0.0);
		 population.get(5).setObjective(1, 1.0);
		 population.get(5).setObjective(2, 0.0);
		 population.get(6).setObjective(0, 1.0);
		 population.get(6).setObjective(1, 0.0);
		 population.get(6).setObjective(2, 0.0);
		 wfghvCalculator2 wfg = new wfghvCalculator2(population,fun);
		 hv = wfg.calculatewfghv();
		 System.out.println("hv= "+ hv);
		
	}

}
