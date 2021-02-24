package jmetal.metaheuristics.moea_c;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
/*
 * ���÷ǵݹ�ķ�����ʵ�ֲ�ξ��࣬���Ľ�����ʱ�临�Ӷ�
 */
public class HierarchicalClusteringInUnitH {
	List<SolutionSet> list = new <SolutionSet>ArrayList();
	private int type = 0;
	public HierarchicalClusteringInUnitH(List list, int type){
		this.list = list;
		this.type = type;
	}
	public List<SolutionSet> clusteringAnalysis(int clusteringSize){
		int size = list.size();
		double minDistance = Double.MAX_VALUE;
		double distance_ = 0.0;
		double[] minDistances = new double[size];//ÿ������֮�����֮��ľ���
		int[] minIndexs = new int[size];//ÿ������֮���������
		int[] index = new int[4];//��ǰ�����Ƶ�������
		index[0] = index[1] = index[2] = index[3] = -1;
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
			SolutionSet sols = (list.get(index[2]).union(list.get(index[3])));
			list.get(index[2]).setRemove(true);	
			list.remove(index[3]);
			list.add(index[3], sols);
			 /*
		     * ���°�index[2]���嵱���Ƕ��������i���������indexs[i];
		     */
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[2] && !list.get(i).isRemove()){
					double min = Double.MAX_VALUE;
		    		int sb = -1;
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j){
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
		     * ���°�index[3]���嵱���Ƕ��������i���������indexs[i];
		     */
			for(int i=0;i<list.size();i++){
				if(minIndexs[i]==index[3] && !list.get(i).isRemove()){
					double min = Double.MAX_VALUE;
		    		int sb = -1;
		    		double ss = 0.0;
		    		for(int j=0;j<list.size();j++){
		    			if(!list.get(j).isRemove() && i!=j){
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
		     * ���µ�ǰ����Ƕȵ����������index;
		     */
			double sAngle = Double.MAX_VALUE;
			for(int p=0;p<list.size();p++){
				if(!list.get(p).isRemove()){
					if(sAngle > minDistances[p]){
						sAngle = minDistances[p];
						index[2] = p;
						index[3] = minIndexs[p]; 	
					}
				}
			}
	
			size--;
		}//while
		/*int sss = 0;
		for(int i=0;i<list.size();i++){
			if(!list.get(i).isRemove()){
				sss++;
			}
		}*/
		/*if(sss != 16){
			System.out.println("sss = " + sss);
			System.exit(0);
		}*/
		Iterator<SolutionSet> iterator = list.iterator();
		while(iterator.hasNext()){
			if(iterator.next().isRemove()){
				iterator.remove();
			}
		}
		return this.list;
	}
	
	/*
     * ����������֮��ĽǶ�ֵ
     */
	public double computeDistance(Solution s1, Solution s2){
		double dis = 0.0;//�������������ĽǶ�
		double p = 2.0;
		for(int i=0; i<list.get(0).get(0).getNumberOfObjectives(); i++){
			dis += Math.pow(Math.abs(s1.getUnitHyperplaneObjective(i) - s2.getUnitHyperplaneObjective(i)), p);
		}
	
	    dis = Math.pow(dis, 1.0/p);
		return dis;
	}//computeDistance

}
