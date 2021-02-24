package jmetal.metaheuristics.moeadds;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
/*
 * ���÷ǵݹ�ķ�����ʵ�ֲ�ξ��࣬���Ľ�����ʱ�临�Ӷ�
 */
public class DividCluster {
	List<SolutionSet> list = new <SolutionSet>ArrayList();
	public DividCluster(List list){
		this.list = list;
	}
	public List<SolutionSet> clusteringAnalysis(int clusteringSize){
		int size = list.size();
		double minAngle = Double.MAX_VALUE;
		double angle = 0.0;
		double[] minAngles = new double[size];//ÿ������֮�����֮��ĽǶ�
		int[] minIndexs = new int[size];//ÿ������֮�������
		int[] index = new int[4];//��ǰ�����Ƶ�������
		index[0] = index[1] = index[2] = index[3] = -1;
		for(int i=0; i<size; i++){
			double min = Double.MAX_VALUE;
			for(int j=0;j<size;j++){
				if(i != j){
					angle = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
					if(min > angle){
						min = angle;
						index[0] = i;
						index[1] = j;
					}
				}
			}
			minAngles[i] =  min;
			minIndexs[i] = index[1];
			if(minAngle > min){
				minAngle = min;
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
		    				ss = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minAngles[i] = min;
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
		    				ss = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
		    				if(min > ss){
			    				min = ss;
			    				sb = j;
			    			}//if
		    			}//if
		    		}//for
		    		minAngles[i] = min;
			    	minIndexs[i] = sb;
				}//if
			}//for
			/*
		     * ���µ�ǰ����Ƕȵ����������index;
		     */
			double sAngle = Double.MAX_VALUE;
			for(int p=0;p<list.size();p++){
				if(!list.get(p).isRemove()){
					if(sAngle > minAngles[p]){
						sAngle = minAngles[p];
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
	public double computeAngle(Solution s1, Solution s2){
		double angle = 0.0;//�������������ĽǶ�
		double distanceToidealPoint1 = s1.getDistanceToIdealPoint();//S1�������ľ���
		double distanceToidealPoint2 = s2.getDistanceToIdealPoint();//S2�������ľ���
		double innerProduc = 0.0; //�����������ڻ�
		for(int i=0; i<list.get(0).get(0).getNumberOfObjectives(); i++){
			innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
		}
		double value = Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2));
		if(value > 1.0){
			//System.out.println(value);
			value = 1.0;
		}
		angle = Math.acos(value);
		return angle;
	}//computeAngle

}
