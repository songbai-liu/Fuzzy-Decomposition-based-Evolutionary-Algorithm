package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.*;

public class wfghvCalculator3 {
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	double[][] pf_ = null;
	double [][] pfMatrix_ = null;
	int fun;
	int number;
	/*public wfghvCalculator3(SolutionSet paretoFront) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;  
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor 
*/	
	/*public wfghvCalculator3(SolutionSet paretoFront,int fun) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;
	    this.fun = fun;
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor
*/	
	public wfghvCalculator3(double[][] paretoFront,int fun, int number) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;
	    this.fun = fun;
	    this.number = number;
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor
	public static double hv2point(double[] Point1, Solution ref){
		  double x=ref.getObjective(0)-Point1[0];
			for (int j=1;j<Point1.length;j++){
				 x = x*(ref.getObjective(j)-Point1[j]);
			}
			return x;
		}
	public double calculatewfghv() throws IOException{
		double[][] sb = new double[pf_.length][number];
		for(int ss=0;ss<pf_.length;ss++){
			for(int j=0;j<number;j++){
				sb[ss][j] = pf_[ss][j];
			}
		}
		
	    double hv; 
	    Solution referencePoint1 = new Solution( number);
		
	    for (int j=0;j<number;j++){
			if(fun == 6){//DTLZ1
				referencePoint1.setObjective(j,0.5);
			}else if(fun<=9){//DTLZ2-DTLZ4
				referencePoint1.setObjective(j,1.0);
			}else if(fun>9&&fun<=11){//DTLZ5 and DTLZ6
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.sqrt(2)/2, j));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun == 12){//MaF7: DTLZ7
				if(j!=number-1){
					referencePoint1.setObjective(j,1.0);
				}else{
					referencePoint1.setObjective(j,2.0*(j+1));
				}
			}else if(fun>12&&fun<=21){//WFG1-WFG9
				referencePoint1.setObjective(j,2.0*(j+1));
			}else if(fun == 22){//MaF1: Inverted DTLZ1
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 23){//MaF2: Concave DTLZ2-BZ
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.cos(Math.PI/8), j+1));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun == 24){//MaF3: Convex DTLZ3
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 25){//MaF4: Inverted badly-scaled DTLZ3
				referencePoint1.setObjective(j,Math.pow(2.0, j+1));
			}else if(fun == 26){//MaF5: Concave-badly scaled DTLZ4
				referencePoint1.setObjective(number-1-j,Math.pow(2.0, j+1));
			}else if(fun == 27){//MaF6: DTLZ5
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.cos(Math.PI/4), j));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun == 28){//MaF7: DTLZ7
				if(j!=number-1){
					referencePoint1.setObjective(j,1.0);
				}else{
					referencePoint1.setObjective(j,2.0*(j+1));
				}
			}else if(fun == 29){//MaF8: Multi-Point Distance Minimization Problem
				double[][] point = new double[number][2];
				point[0][0] = 0.0;
				point[0][1] = 1.0;	
				double arc = 2*Math.PI/number;
				for (int i = 1; i < number; i++){
					point[i][0] = point[0][0] - Math.sin(arc*i);
					point[i][1] = point[0][1] - 1.0 + Math.cos(arc*i);
				}
				double maxValue = Double.MIN_VALUE;
				for(int i=0;i<number;i++){
					for(int s=i+1;s<number;s++){
						double value = Math.pow(point[i][0]-point[j][0], 2) + Math.pow(point[i][1]-point[j][1], 2);
						if(value > maxValue){
							maxValue = value;
						}
					}	
				}
				referencePoint1.setObjective(j,Math.sqrt(maxValue));
			}else if(fun == 30){//MaF9: Multi-Line Distance Minimization Problem
				double[][] point = new double[number][2];
				point[0][0] = 0.0;
				point[0][1] = 1.0;	
				double arc = 2*Math.PI/number;
				for (int i = 1; i < number; i++){
					point[i][0] = point[0][0] - Math.sin(arc*i);
					point[i][1] = point[0][1] - 1.0 + Math.cos(arc*i);
				}
				
				double[] k = new double[number];
				double[] f = new double[number];
				for (int i = 0; i < number-1; i++){
					k[i] = (point[i+1][1]-point[i][1])/(point[i+1][0]-point[i][0]);
					f[i] = Math.abs(point[0][1]-k[i]*point[0][0]+k[i]*point[i][0]-point[i][1])/Math.sqrt(1+Math.pow(k[i], 2));
				}
				k[number-1]=(point[number-1][1]-point[0][1])/(point[number-1][0]-point[0][0]);
				f[number-1] = Math.abs(point[0][1]-k[number-1]*point[0][0]+k[number-1]*point[0][0]-point[0][1])
						/Math.sqrt(1+Math.pow(k[number-1], 2));
				
				double maxValue = Double.MIN_VALUE;
				for(int i=0;i<number;i++){
					if(f[i] > maxValue){
						maxValue = f[i];
					}	
				}
				referencePoint1.setObjective(j,Math.sqrt(maxValue));
			}else if(fun>=31&&fun<=33){//MaF10-MaF12: WFG1, WFG2, and WFG9
				referencePoint1.setObjective(j,2.0*(j+1));
			}else if(fun == 34){//MaF13: PF7
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 35){//mDTLZ1
			    referencePoint1.setObjective(j,0.5);
			}else if(fun>=36&&fun<=38){//mDTLZ2-mDTLZ4
				referencePoint1.setObjective(j,1.0);
			}else if(fun == 39){//
				
			}
		}
		
		//NORMALIZATION
		for (int j=0;j<sb.length;j++)
			for(int k=0;k<number;k++)
				sb[j][k] = (sb[j][k])/(1.1*referencePoint1.getObjective(k));
		        //sb.get(j).setObjective(k, sb.get(j).getObjective(k)/(referencePoint1.getObjective(k)) );
		//SolutionSet invertedFront;
		//invertedFront = utils_.invertedFront(pf_,number);
		for (int j=0;j<sb.length;j++)
			for(int k=0;k<number;k++)
				if(sb[j][k]>1.0){
					for(int s=0;s<number;s++){
						sb[j][s]  = 1.0;
					}
					break;
				}
		for (int j=0;j<number;j++)
		    referencePoint1.setObjective(j,1.0);
		
		if (sb.length == 0){
			hv = 0.0;
		}else if(sb.length == 1){
			hv = hv2point(sb[0], referencePoint1);
		}else if (sb.length == 2){
			double [] mid = new double [number];
			for (int j=0;j<number;j++){
				mid[j] = Math.max(sb[0][j], sb[1][j]);
			}
			double[] midp = new double[number];
			for(int i=0;i<number;i++){
				midp[i] = mid[i];
			}
			hv = hv2point(sb[0],referencePoint1)+hv2point(sb[1],referencePoint1)-hv2point(midp,referencePoint1);
			
		}else if (sb.length==3){
			double [] w01 = new double [number];
			double [] w02 = new double [number];
			double [] w12 = new double [number];
			double [] w012 = new double [number];
			for (int j=0;j<number;j++){
				w01[j] = Math.max(sb[0][j], sb[1][j]);
				w02[j] = Math.max(sb[0][j], sb[2][j]);
				w12[j] = Math.max(sb[1][j], sb[2][j]);
			}
			for (int j=0;j<number;j++){
				w012[j] = Math.max(w02[j], sb[1][j]);
			}
			hv = hv2point(sb[0],referencePoint1)+hv2point(sb[1],referencePoint1)+hv2point(sb[2],referencePoint1)
					-hv2point(w01,referencePoint1)-hv2point(w02,referencePoint1)-hv2point(w12,referencePoint1)+hv2point(w012,referencePoint1);
		}else{
			WFGHV1 wfghv = new WFGHV1(number, sb.length) ;
			Front1 front = new Front1(sb.length, number, sb);
		    hv = wfghv.getHV(front,referencePoint1);
		}
		return hv;
  } // CalculateHypervolume
}
