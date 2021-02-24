package jmetal.metaheuristics.moea_c;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;
/*
 * 采用非递归的方法来实现层次聚类，大大的降低了时间复杂度
 */
public class HierarchicalClusteringForVectors {
	List<SolutionSet> list = new <SolutionSet>ArrayList();
	private int type = 6;
	
	private double[][] vectors;
	
	public HierarchicalClusteringForVectors(List list){
		this.list = list;
	}
	
	public HierarchicalClusteringForVectors(List list, int type){
		this.list = list;
		this.type = type;
	}
	public List<SolutionSet> clusteringAnalysis(int clusteringSize){
		int size = list.size();
		double minDistance = Double.MAX_VALUE;
		double distance_ = 0.0;
		double[] minDistances = new double[size];//每个类与之最近类之间的距离
		int[] minIndexs = new int[size];//每个类与之最近的类标号
		int[] index = new int[4];//当前最相似的两个类
		index[0] = index[1] = index[2] = index[3] = -1;
		boolean[] isFull = new boolean[size]; 
		for(int i=0; i<size; i++){
			double min = Double.MAX_VALUE;
			for(int j=0;j<size;j++){
				if(i != j){
					distance_ = computeDistance(list.get(i).getCentroid(),list.get(j).getCentroid());
					if(min > distance_){
						min = distance_;
						index[0] = i;
						index[1] = j;
					}
				}
			}
			minDistances[i] =  min;
			minIndexs[i] = index[1];
			if(minDistance > min){
				minDistance = min;
				index[2] = index[0];
				index[3] = index[1];
			}
		}
		while(size > clusteringSize){
			if(index[2] == -1 || index[3] == -1){
				System.out.println("index[2] = " +index[2]+" ,index[3] = "+index[3]);
				System.exit(0);
			}
			
			
			SolutionSet sols = (list.get(index[2]).union(list.get(index[3])));
			list.get(index[2]).setRemove(true);	
			list.remove(index[3]);
			list.add(index[3], sols);
			if(sols.size() >= type){
				isFull[index[3]]= true;
			}
			 /*
		     * 更新把index[2]个体当作角度最近个体i的最近个体indexs[i];
		     */
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[2] && !list.get(i).isRemove() && !isFull[i]){
					double min = Double.MAX_VALUE;
		    		int sb = -1;
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j && !isFull[j]){
		    				ss = computeDistance(list.get(i).getCentroid(),list.get(j).getCentroid());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minDistances[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			 /*
		     * 更新把index[3]个体当作角度最近个体i的最近个体indexs[i];
		     */
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[3] && !list.get(i).isRemove() && !isFull[i]){
					double min = Double.MAX_VALUE;
		    		int sb = -1;
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j && !isFull[j]){
		    				ss = computeDistance(list.get(i).getCentroid(),list.get(j).getCentroid());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minDistances[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			/*
		     * 更新当前最近角度的两个个体的index;
		     */
			double sAngle = Double.MAX_VALUE;
			for(int p=0;p<list.size();p++){
				if(!list.get(p).isRemove() && !isFull[p]){
					if(sAngle > minDistances[p]){
						sAngle = minDistances[p];
						index[2] = p;
						index[3] = minIndexs[p]; 	
					}
				}
			}
			size--;
			//size--;
		}//while
		
		Iterator<SolutionSet> iterator = list.iterator();
		while(iterator.hasNext()){
			if(iterator.next().isRemove()){
				iterator.remove();
			}
		}
		return this.list;
	}
	
	public double[][] getVectors(int clusterSize){
		List<SolutionSet> clusters = this.clusteringAnalysis(clusterSize);
		int number = list.get(0).get(0).getNumberOfObjectives();
		this.vectors = new double[clusterSize][number];
		
		for(int k=0; k<number;k++){
    		double minClustering2Axis = 1.0e+30;
      		int minClustering2AxisID = -1;
      		for(int i=0;i<list.size();i++){
      			SolutionSet sols = list.get(i);
      			if(sols.size() == 0){
      				System.out.println("SolsSize1 = "+sols.size());
      				System.exit(0);
      			}
      			
      			double angle1 = Math.acos(Math.abs(sols.getCentroid().getNormalizedObjective(k)/sols.getCentroid().getDistanceToIdealPoint()));
      			//System.out.println(angle1);
      			if(angle1 < minClustering2Axis){
      				minClustering2Axis = angle1;
      				minClustering2AxisID = i;
      			}//if
      		}//for
      		if(minClustering2AxisID == -1){
  				System.out.println("minClustering2AxisID = -1");
  				System.exit(0);
  			}
      		double minSolution2Axis = 1.0e+30;
      		int minSolution2AxisID = -1;
      		for(int j=0;j<list.get(minClustering2AxisID).size();j++){
      			Solution sol = list.get(minClustering2AxisID).get(j);
      			double ang = Math.acos(list.get(minClustering2AxisID).get(j).getNormalizedObjective(k)/list.get(minClustering2AxisID).get(j).getDistanceToIdealPoint());
      			if(ang < minSolution2Axis){
      				minSolution2Axis = ang;
      				minSolution2AxisID = j;
      			}
      		}//for
      		for(int i=0;i<number;i++){
				vectors[k][i] = list.get(minClustering2AxisID).get(minSolution2AxisID).getNormalizedObjective(i);
			}
            list.remove(minClustering2AxisID);
		}
		Iterator<SolutionSet> iterator = clusters.iterator();
		int p = number;
		while(iterator.hasNext()){
			Solution sb = iterator.next().getCentroid();
			for(int i=0;i<number;i++){
				vectors[p][i] = sb.getNormalizedObjective(i);
			}
			p++;
		}
		return this.vectors;
	}
	
	/*
     * 求两个个体之间的角度值
     */
	public double computeDistance(Solution s1, Solution s2){
		double dis = 0.0;
		double p = 2.0;
		for(int i=0; i<list.get(0).get(0).getNumberOfObjectives(); i++){
			dis += Math.pow(Math.abs(s1.getUnitHyperplaneObjective(i) - s2.getUnitHyperplaneObjective(i)), p);
		}
	
	    dis = Math.pow(dis, 1.0/p);
		return dis;
	}//computeDistance

}
