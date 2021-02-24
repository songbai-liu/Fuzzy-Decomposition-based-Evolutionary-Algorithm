package jmetal.metaheuristics.maoeac;

import java.util.ArrayList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
/*
 * ���õݹ�ķ�����ʵ�ֲ�ξ��࣬ʱ�临�Ӷ��ر��
 */
public class HierarchicalClustering {
	List<SolutionSet> list = new <SolutionSet>ArrayList();
	public HierarchicalClustering(List list){
		this.list = list;
	}
	public List<SolutionSet> clusteringAnalysis(int clusteringSize){
		int size = list.size();
		double minAngle = Double.MAX_VALUE;
		double angle = 0.0;
		double[] minAngles = new double[size];//ÿ������֮�����֮��ĽǶ�
		int[] minIndexs = new int[size];//ÿ������֮�������
		int[] index = new int[2];//��ǰ�����Ƶ�������
		index[0] = index[1] = -1;
		for(int i=0; i<size; i++){
			for(int j=0;j<size;j++){
				if(i != j){
					angle = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
					if(minAngle > angle){
						   minAngle = angle;
						   index[0] = i;
						   index[1] = j;
					}
				}
			}
			minAngles[i] =  minAngle;
			minIndexs[i] = index[1];
		}
		if(size > clusteringSize){
			SolutionSet sols = (list.get(index[0]).union(list.get(index[1])));
			if(index[0] > index[1]){
				list.remove(index[0]);
				list.remove(index[1]);
			}else{
				list.remove(index[1]);
				list.remove(index[0]);
			}
			list.add(sols);
		}else{
			return this.list;
		}
		clusteringAnalysis(clusteringSize);
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
		angle = Math.acos(Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2)));
		return angle;
	}//computeAngle

}
