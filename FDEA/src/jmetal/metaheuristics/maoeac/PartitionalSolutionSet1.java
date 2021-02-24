package jmetal.metaheuristics.maoeac;

import jmetal.core.SolutionSet;

public class PartitionalSolutionSet1 {
	SolutionSet solutionSet_;
	int numberOfClusters;
	int[][] neighbor;
	
	public PartitionalSolutionSet1(SolutionSet solutionSet_,int numberOfObjectives){
		this.solutionSet_ = solutionSet_;
		this.numberOfClusters = numberOfObjectives;
		neighbor = new int[numberOfObjectives][solutionSet_.size()/(numberOfObjectives+1)];
	}
	
	public SolutionSet[] partitional(){
		SolutionSet[] sols = new SolutionSet[numberOfClusters+1];
		int sb = solutionSet_.size() / (numberOfClusters+1);
		sols[numberOfClusters] = new SolutionSet(sb);
		SolutionSet solset = new SolutionSet(solutionSet_.size());
		for(int p=0;p<solutionSet_.size(); p++){
			solset.add(solutionSet_.get(p));
		}
		double[] ux = new double[solset.size()];
		int[] idux = new int[solset.size()];
		for(int q=0;q<solset.size();q++){
			ux[q] = Math.acos(Math.abs(solset.get(q).getSumValue()/(solset.get(q).getDistanceToIdealPoint()*Math.sqrt(solset.get(0).getNumberOfObjectives()))));
		    idux[q] = q;
		}
		Utils.minFastSort(ux, idux, solset.size(), sb);
		boolean[] flag = new boolean[solset.size()];
		for(int k=0;k<sb;k++){
			flag[idux[k]] = true;
		}
		SolutionSet solsets = new SolutionSet(solset.size()-sb);
		for(int s=0;s<solset.size();s++){
			if(flag[s] == true){
				sols[numberOfClusters].add(solset.get(s));
			}else{
				solsets.add(solset.get(s));
			}
		}
		
		int sbSize = solsets.size() / numberOfClusters;
		boolean[] isAsocciated = new boolean[solsets.size()];
		for(int i=0;i<numberOfClusters;i++){
			sols[i] = new SolutionSet(sbSize);
		}
		int[] permutation = new int[numberOfClusters];
		Utils.randomPermutation(permutation,numberOfClusters);
		for(int i=0;i<numberOfClusters;i++){
			int n = permutation[i];
			int currentSize = solsets.size();
			for(int p=0;p<solsets.size();p++){
				if(isAsocciated[p]){
					currentSize--;
				}
			}
			double[] x = new double[currentSize];
			int[] idx = new int[currentSize];
			int t = 0;
			for(int j=0;j<solsets.size();j++){
				if(!isAsocciated[j]){
					x[t] = Math.acos(Math.abs(solsets.get(j).getNormalizedObjective(n)/solsets.get(j).getDistanceToIdealPoint()));
					idx[t] = j;
					t++;
				}
			}
			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, currentSize, sbSize);

			System.arraycopy(idx, 0, neighbor[n], 0, sbSize);
			for(int k=0;k<neighbor[n].length;k++){
				isAsocciated[neighbor[n][k]] = true;
			}
		}
		for(int i=0;i<numberOfClusters;i++){
			for(int j=0;j<sbSize;j++){
				sols[i].add(solsets.get(neighbor[i][j]));
			}
		}
		return sols;
	}

}
