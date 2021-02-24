package jmetal.metaheuristics.maoeac;

import java.util.ArrayList;
import java.util.List;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.comparators.GaussianLocalDensityComparator;

public class PartitionalClustering {
	private double dtheata;
	SolutionSet solutionSet;
	public PartitionalClustering(SolutionSet sols){
		this.solutionSet = sols;
		this.dtheata = Math.PI*8/(sols.size());
	}
	
	public PartitionalClustering(SolutionSet sols,double dtheata){
		this.solutionSet = sols;
		this.dtheata = dtheata;
	}
	
	 /*
     * 求两个个体之间的角度值
     */
	public double computeAngle(Solution s1, Solution s2){
		double angle = 0.0;//所求两个向量的角度
		double distanceToidealPoint1 = s1.getDistanceToIdealPoint();//S1到理想点的距离
		double distanceToidealPoint2 = s2.getDistanceToIdealPoint();//S2到理想点的距离
		double innerProduc = 0.0; //两个向量的内积
		for(int i=0; i<solutionSet.get(0).getNumberOfObjectives(); i++){
			innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
		}
		angle = Math.acos(Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2)));
		return angle;
	}//computeAngle
	
	public void computeLocalDensity(SolutionSet union){
		double[] ld = new double[union.size()];
		
		for(int i=0;i<union.size();i++){
			ld[i] = 0;
			Solution s1 = union.get(i);
			for(int j=0;j<union.size();j++){
				Solution s2 = union.get(j);
				if(i != j){
					double angle = computeAngle(s1,s2);
					if (angle <= dtheata){
						ld[i] += Math.pow(Math.E, -(Math.pow(angle/dtheata, 2)));	
						
					}					
				}
			}
			//System.out.print(ld[i]+" ");
		}
		for(int s=0;s<union.size();s++){
			union.get(s).setGaussianLocalDensity(ld[s]);
		}
	}
	
	public void computeCenterDistance(SolutionSet union){
		double[] cd = new double[union.size()];
		union.sort(new GaussianLocalDensityComparator());
		for(int i=1;i<union.size();i++){
			double minAng = 1.0e+30;
			Solution s1 = union.get(i);
			for(int j=0;j<i;j++){
				Solution s2 = union.get(j);
				double angle = computeAngle(s1,s2);
				if (angle < minAng){
					minAng = angle;					
				}					
			}
			cd[i] = minAng;
		}
		double maxAng = Double.MIN_VALUE;
		for(int k=1;k<union.size();k++){
			if(maxAng < cd[k]){
				maxAng = cd[k];
			}
		}
		cd[0] = maxAng;
		
		for(int s=0;s<union.size();s++){
			union.get(s).setCenterDistance(cd[s]);
		}
	}
	
	public List<SolutionSet> clusteringAnalysis(){
		computeLocalDensity(solutionSet);
		computeCenterDistance(solutionSet);
		double[] p = new double[solutionSet.size()];
		double[] q = new double[solutionSet.size()];
		double[] r = new double[solutionSet.size()];
		
		double sump = 0.0;
	    double avep = 0.0;
	    double sumq = 0.0;
	    double aveq = 0.0;
		
		for(int i=0;i<solutionSet.size();i++){
	    	p[i] = solutionSet.get(i).getGaussianLocalDensity();
	    	q[i] = solutionSet.get(i).getCenterDistance();
	    	sump += p[i];
	    	sumq += q[i];
	    	r[i] = p[i]*q[i];
	    	solutionSet.get(i).setCrowdingDistance(r[i]);
	    }
		
		avep = sump/solutionSet.size();
		aveq = sumq/solutionSet.size();
		if(aveq < dtheata){
			System.out.println("average center distance is small d = "+ aveq + " detheata = "+dtheata);
		    //aveq = dtheata;
		}
		  
		SolutionSet sols = new SolutionSet();
		boolean[] isAdded = new boolean[solutionSet.size()];
		for(int j=0;j<solutionSet.size();j++){
		    if( (p[j] > 1.0*avep && q[j] > 1.0*aveq)){
		    	sols.add(solutionSet.get(j));
		    	isAdded[j] = true;
		    }
		}
		int[] ids = new int[solutionSet.size()];
		    
		//Associate All Solutions with a Clustering Center
		for(int k=0;k<solutionSet.size();k++){
		    if(!isAdded[k]){//is not the clustering center
		    	Solution s1 = solutionSet.get(k);
		    	int minInd = -1;
		    	double minAngle = 1.0e+30;
		    	for(int l=0;l<solutionSet.size();l++){
		    		if(isAdded[l]){//is the clustering center
		    			Solution s2 = solutionSet.get(l);
			    		double angle = computeAngle(s1,s2);
			    		if (angle < minAngle){
							minAngle = angle;
							minInd = l;
						}
		    		}
		    		if(minAngle < dtheata){
			    		ids[k] = minInd;
			    	}else{
			    		ids[k] = -1;
			    	}
			    }
			 }else{
			    ids[k] = k;
			 }
		}
		//solutionSet.sort(new CrowdingDistanceComparator());
		    
	    boolean[] isRemoved = new boolean[solutionSet.size()];
		SolutionSet[] clusterings = new SolutionSet[sols.size()];
		for(int ss=0;ss<sols.size();ss++){
		    clusterings[ss] = new SolutionSet();
		}
		int c=0;
		    //int maxT = -1;
		for(int n=0;n<solutionSet.size();n++){
		    if(isAdded[n]){
		    	for(int m=0;m<solutionSet.size();m++){
			    	if(ids[m] == n){
			    		clusterings[c].add(solutionSet.get(m));
			    		isRemoved[m] = true;
			    	}
			    }
		    	c++;
		    }
		}
		
		List<SolutionSet> list = new <SolutionSet>ArrayList();
		
		for(int i=0;i<sols.size();i++){
			list.add(clusterings[i]);
		}
		
		for(int j=0;j<solutionSet.size();j++){
			if(!isRemoved[j]){
				SolutionSet ss = new SolutionSet();
				ss.add(solutionSet.get(j));
				list.add(ss);
			}
		}
		
		return list;
	}
	
	

}
