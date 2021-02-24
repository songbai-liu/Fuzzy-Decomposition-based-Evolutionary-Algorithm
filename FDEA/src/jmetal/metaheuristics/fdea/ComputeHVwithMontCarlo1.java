package jmetal.metaheuristics.fdea;
import jmetal.util.RandomGenerator;

public class ComputeHVwithMontCarlo1 {
	
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	 /**
	 * Constructor
	 * Creates a new instance of MultiDelta
	 */
	public ComputeHVwithMontCarlo1(){
		 utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	}
	
	 /***
	   * estimate the hypervolume 
	   * */
       public double estimateHypervolumeMonteCarlo(double p,
		  double[] maximumValues,  double[] minimumValues,  
		 long samplingNumberM, int numberOfObjectives) throws ClassNotFoundException {
		/***
		 * estimate hypervolume con
		 * record the value in the "fitness" Field of each solution
		 */
	   
		int nObj = numberOfObjectives; 
		double sp = p;

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
			double sum1 = 0;
			double sum2 = 0;
			for(int i = 0; i < nObj; i ++){
				sum1 += Math.pow(s[i], sp);
				sum2 += s[i];
			}
			if(sp == 1.0){
				break;
			}else if(sp > 1.0){
				if(sum1 <= 1.0 && sum2 >= 1.0){
					k++;
				}
			}else{
				if(sum1 >= 1.0 && sum2 <= 1.0){
					k++;
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
