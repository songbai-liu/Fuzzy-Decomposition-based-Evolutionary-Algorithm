package jmetal.metaheuristics.fdea;
import jmetal.util.RandomGenerator;

public class ComputeHVwithMontCarlo {
	
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	 /**
	 * Constructor
	 * Creates a new instance of MultiDelta
	 */
	public ComputeHVwithMontCarlo(){
		 utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	}
	
	/* 
	   returns true if 'point1' dominates 'points2' with respect to the 
	   to the first 'noObjectives' objectives 
	   */
	 public boolean  dominates(double  point1[], double  point2[], int  noObjectives) {
	    int  i;
	    int  betterInAnyObjective;

	    betterInAnyObjective = 0;
	    for (i = 0; i < noObjectives && point1[i] <= point2[i]; i++)
	      if (point1[i] < point2[i]) 
	      	betterInAnyObjective = 1;
	    
	    return ((i >= noObjectives) && (betterInAnyObjective>0));
	  } //Dominates
	
	 /***
	   * estimate the hypervolume 
	   * */
       public double estimateHypervolumeMonteCarlo(double[][] solutionset,
		  double[] maximumValues,  double[] minimumValues,  
		 long samplingNumberM, int numberOfObjectives) throws ClassNotFoundException {
		/***
		 * estimate hypervolume con
		 * record the value in the "fitness" Field of each solution
		 */
	   
		int nObj = numberOfObjectives;  

     //		private void refreshReferenceSet(SolutionSet rSet, SolutionSet population) 
		/*determine sampling box S*/
		SamplingBox S = new SamplingBox(nObj);
		
		S.setBoundsValue(maximumValues, minimumValues);
		S.setVValues();
		double V = S.getV();
		//System.out.println("V==" + V);
		/*reset fitness assignment*/
		/*perform sampling*/
		double k = 0;
		for(int j = 0 ;j < samplingNumberM; j ++ ) {
			double[] s = S.getOneRandomSample();
			for(int i = 0; i < solutionset.length; i ++){
				double[] point1 = solutionset[i];
				if(dominates(point1, s, numberOfObjectives)){
					k ++;
					break;
				}
			}
			 
		}//for sampling	 
		//System.out.println("k= "+k);
		double hvOfSolutionSet = V*(double)k/(double)samplingNumberM;
		return hvOfSolutionSet;
	}


	private class SamplingBox{
		double[] lowerBounds;
		double[] upperBounds;
		int nObj;
		private RandomGenerator myRandom = new RandomGenerator(0); 
		double V = 1;
		public SamplingBox(int nObj){
			this.nObj = nObj;
			this.lowerBounds = new double[nObj];
			this.upperBounds = new double[nObj];
		}
		public double[] getOneRandomSample(){
			double[] ss = new double[this.nObj];
			for(int i = 0; i < this.nObj; i ++){
				double length = this.upperBounds[i] - this.lowerBounds[i];
				ss[i] =  this.myRandom.nextDouble() *  length + this.lowerBounds[i];
			}
			return ss;
		}
		public void setBoundsValue(double[] maximumValues, double[] minimumValues){
			for(int i = 0; i< nObj; i ++) { 
				lowerBounds[i] = minimumValues[i] ;
				upperBounds[i] = maximumValues[i] ;
			}
			
		}
		public void setVValues(){
			for(int i = 0; i < nObj; i ++) {
				double length = this.upperBounds[i] - this.lowerBounds[i];
				//double length = 1.0 - this.lowerBounds[i];
				if(length > 0){
					this.V *=length;
				}
				else {
					this.V = 0;
					//this.V *= -length;
				}
			}
		}
		public double getV(){
			return V;
		}
			
	}

}
