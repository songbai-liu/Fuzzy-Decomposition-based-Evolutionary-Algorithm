package jmetal.metaheuristics.fdea;

import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import jmetal.core.Problem;
import jmetal.problems.ProblemFactory;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

public class PSample_main{
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *             paretoFrontFile
	 */
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区    
	      
	      NumberFormat nf = NumberFormat.getInstance();
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(nf.format(GD[i])+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void printGD(String path,double[] GD, double[] HV){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区         
	      NumberFormat nf = NumberFormat.getInstance();
	      DecimalFormat df = new DecimalFormat("0.00000000");
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(nf.format(GD[i])+" ");//写到缓冲区
	        bw.write(df.format(HV[i])+" ");
	    	//bw.write(GD[i]+" ");//写到缓冲区
		    //bw.write(HV[i]+" ");
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	public static void main(String args[]) throws JMException, ClassNotFoundException, SecurityException, IOException{
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAIII_main.log");
		logger_.addHandler(fileHandler_);
		
		for(int fun=3;fun<=4;fun++){
		  double[] p_samples;
		  double[] HV_Samples;
		  
		  Problem problem = null; // The problem to solve
		  QualityIndicator indicators;// Object to get quality indicators
			indicators = null;
			
			if (args.length == 1) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
				indicators = new QualityIndicator(problem, args[1]);
			} // if
			else { // Default problem
			    if(fun==1){
				      problem = new DTLZ1("Real",7,2);
			    }
			    if(fun==2){
				      problem = new DTLZ1("Real",7,3);
			    }
			    if(fun==3){
				      problem = new DTLZ1("Real",7,5);
			    }
			    if(fun==4){
				      problem = new DTLZ1("Real",7,8);
			    }
			    if(fun==5){
				      problem = new DTLZ1("Real",7,10);
			    }
			    if(fun==6){
				      problem = new DTLZ1("Real",7,15);
			    }
			}
			
			int numb = problem.getNumberOfObjectives();
			for(int k=1;k<=numb;k++){
				int s = 1;
				for(int t=0;t<k;t++){
					s*=(numb-t);
				}
			}
	    	int K1 = 9;
	    	int K2 = 20;
	    	int t = K1+K2;
	    	p_samples = new double[t];
	    	
	    	for(int k=0; k<K1; k++){
	    		p_samples[k] = 1 - 0.1*k;
	    	}
	    	for(int k=0; k<K2; k++){
	    		p_samples[K1+k] = 1 + 0.2*(k+1);
	    	}
	    	
	    	double[] referencePoint = new double[numb];
			double[] referencePoint1 = new double[numb];
			
			for(int i=0;i<numb;i++){
				referencePoint1[i] = 0.0;
				referencePoint[i] = 1.0;
			}
			long samplingNumberM = 50000000;
			double S = Math.sqrt(numb);
			/*for(int m=2;m<numb;m++){
				S = S/m;
			}*/
			S = S/(numb-1);
	    	HV_Samples = new double[t];
	    	
	    	ComputeHVwithMontCarlo1 computeHV = new ComputeHVwithMontCarlo1();
	    	
	    	for(int i=0;i<t;i++){
	    	   HV_Samples[i] = computeHV.estimateHypervolumeMonteCarlo(p_samples[i], referencePoint, referencePoint1, samplingNumberM, numb);
	    	   if(p_samples[i] >= 1.0){
	    		   HV_Samples[i] = HV_Samples[i]/S;
	    	   }else{
	    		   HV_Samples[i] = -HV_Samples[i]/S;
	    	   }
	    	  
	    	}
	    	printGD("PSamples_"+problem.getNumberOfObjectives()+"Obj_"+"_H.txt",p_samples,HV_Samples);
		}
		
		
		
  }//main
}
