package jmetal.metaheuristics.maoeac;

import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.io.*;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.thetadea.ThetaDEA;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.learning.LearningFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.LSMOP.LSMOP1;
import jmetal.problems.LSMOP.LSMOP2;
import jmetal.problems.LSMOP.LSMOP5;
import jmetal.problems.LSMOP.LSMOP9;
import jmetal.problems.MaF.MaF1;
import jmetal.problems.MaF.MaF10;
import jmetal.problems.MaF.MaF11;
import jmetal.problems.MaF.MaF12;
import jmetal.problems.MaF.MaF13;
import jmetal.problems.MaF.MaF2;
import jmetal.problems.MaF.MaF3;
import jmetal.problems.MaF.MaF4;
import jmetal.problems.MaF.MaF5_Concave;
import jmetal.problems.MaF.MaF5_Convex;
import jmetal.problems.MaF.MaF6;
import jmetal.problems.MaF.MaF7;
import jmetal.problems.MaF.MaF8;
import jmetal.problems.MaF.MaF9;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
import jmetal.problems.mDTLZ.mDTLZ1;
import jmetal.problems.mDTLZ.mDTLZ2;
import jmetal.problems.mDTLZ.mDTLZ3;
import jmetal.problems.mDTLZ.mDTLZ4;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.util.Configuration;
import jmetal.util.JMException;

public class MLMOEA_main_10obj {
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
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
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
		
		for(int fun=41;fun<=42;fun++){
		int runtimes=4;
		double[] GDarray = new double[runtimes];
		double[] spreadarray = new double[runtimes];
		double[] IGDarray=new double[runtimes];	
		double[] HVarray = new double[runtimes];
		double[] times = new double[runtimes];
		long Execution_time=0;
		
		Problem problem = null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection; //Selection operator
		Operator learning1; //Learning operator
		Operator learning2; //Learning operator
		
		HashMap parameters; // Operator parameters
		
		Solution referencePoint;
		
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
		  	      problem = new ZDT1("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT1_501.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==2){
		  	      problem = new ZDT2("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT2_501.txt" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==3){
		  	      problem = new ZDT3("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT3_269.txt" ) ;
		  	    	}
		  	if(fun==4){
			      problem = new ZDT4("Real",10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT4_501.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==5){
			      problem = new ZDT6("Real",10);
			      
			      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT6_774.txt" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==6){
				//problem = new DTLZ1("Real",10,2);
			      problem = new DTLZ1("Real",14,10);
			      indicators = new QualityIndicator(problem) ;
			      //indicators = new QualityIndicator(problem,"F:\\sbMaOP\\SbMaOP\\truePF\\DTLZ1_M10.dat" ) ;
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",19,10);
		  		indicators = new QualityIndicator(problem) ;
			      //indicators = new QualityIndicator(problem,"F:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",19,10);
			      indicators = new QualityIndicator(problem) ;
			      //indicators = new QualityIndicator(problem,"F:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",19,10);
			      indicators = new QualityIndicator(problem) ;
			      //indicators = new QualityIndicator(problem,"F:\\sbMaOP\\SbMaOP\\truePF\\DTLZ2_M10.dat" ) ;
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",19,10);
			      indicators = new QualityIndicator(problem) ;
			      //indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ5_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",19,10);
			      indicators = new QualityIndicator(problem) ;
			      //indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\DTLZ5_M10.dat" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",29,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",28,20,15);
		  	      
		  	      indicators = new QualityIndicator(problem) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",28,20,15);
		  	      
		  	      indicators = new QualityIndicator(problem) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",28,20,15);
		  	      
		  	      indicators = new QualityIndicator(problem) ;
		  	    	}
			if(fun==16){
			      problem = new WFG4("Real",28,20,15);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==17){
			      problem = new WFG5("Real",28,20,15);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
		  	if(fun==18){
		  	      problem = new WFG6("Real",28,20,15);
		  	      
		  	      indicators = new QualityIndicator(problem) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
			      problem = new WFG7("Real",28,20,15);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==20){
			      problem = new WFG8("Real",28,20,15);
			      
			      indicators = new QualityIndicator(problem);
			    	}
			if(fun==21){
			      problem = new WFG9("Real",28,20,15);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==22){
			      problem = new MaF1("Real",19,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==23){
			      problem = new MaF2("Real",19,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==24){
			      problem = new MaF3("Real",19,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==25){
			      problem = new MaF4("Real",19,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==26){
			      //problem = new MaF5_Convex("Real",14,5);
			      problem = new MaF5_Concave("Real",19,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==27){
			      problem = new MaF6("Real",19,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==28){
			      problem = new MaF7("Real",29,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==29){
			      problem = new MaF8("Real",2,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==30){
			      problem = new MaF9("Real",2,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==31){
			      problem = new MaF10("Real",18,20,10);//WFG1
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==32){
			      problem = new MaF11("Real",18,20,10);//WFG2
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==33){
			      problem = new MaF12("Real",18,20,10);//WFG9
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==34){
			      problem = new MaF13("Real",5,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			 if(fun==35){
			      problem = new mDTLZ1("Real",20,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==36){
			      problem = new mDTLZ2("Real",20,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==37){
			      problem = new mDTLZ3("Real",20,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==38){
			      problem = new mDTLZ4("Real",20,10);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==39){
			      problem = new LSMOP1("Real",300,5);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP1_5D_20000.txt");
			    	}
			if(fun==40){
			      problem = new LSMOP2("Real",300,10);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP2_10D_30000.txt");
			    	}
			if(fun==41){
			      problem = new LSMOP5("Real",300,10);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP5_2D_5000.txt");
			    	}
			if(fun==42){
			      problem = new LSMOP9("Real",300,10);
			      indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP9_2D_5000.txt");
			    	}
		} // else

		algorithm = new MaOEA_C(problem);
		algorithm.setInputParameter("populationSize", 440);
		algorithm.setInputParameter("maxGenerations",4500000);
		
		/*algorithm.setInputParameter("T",20);
		algorithm.setInputParameter("delta", 0.8);*/
		
		// Mutation and Crossover for Real codification
		parameters = new HashMap();
		parameters.put("probability", 1.0);
		parameters.put("distributionIndex", 30.0);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
				parameters);

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation",
				parameters);

		parameters = null;
		selection = SelectionFactory.getSelectionOperator("RandomSelection",
				parameters);
		
	/*	parameters = new HashMap();
		parameters.put("learningProbability", 1.0);
		parameters.put("learningRate",  0.9);
		learning1 = LearningFactory.getLearningOperator("IntegratedLearning", parameters);
		
		parameters = new HashMap();
		parameters.put("learningProbability", 1.0);
		parameters.put("learningRate",  0.9);
		learning2 = LearningFactory.getLearningOperator("DeepLearning", parameters);*/

		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);
		//algorithm.addOperator("learning1", learning1);
		//algorithm.addOperator("learning", learning2);
		
		for(int i=0;i<runtimes;i++){
		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		Execution_time+=(System.currentTimeMillis() - initTime);
		//times[i] = (double)(System.currentTimeMillis() - initTime)/1000;
		population.printObjectivesToFile("MaOEAC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
		/*if(fun>=29 && fun<=30){
			population.printVariablesToFile("MaOEAC_Variables_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"run_"+ (i+1)  + ".txt");
		}*/
		//IGDarray[i]=indicators.getIGD1(population);
		//wfghvCalculator1 wfg = new wfghvCalculator1(population,fun);
		//HVarray[i] = wfg.calculatewfghv();
		}
		//printGD("MaOEA_HC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_runtimes.txt",times);
		//printGD("MaOEAC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_HV.txt",HVarray);
		//double sumIGD=0;
		double sumHV=0.0;
		for(int i=0;i<runtimes;i++){
			  //sumIGD+=IGDarray[i];
			  sumHV+=HVarray[i];
			}
		logger_.info("Total execution time: " + Execution_time + "ms");
	    //System.out.println("avrIGD-fun"+fun+"= "+sumIGD/runtimes);
		//System.out.println();
		System.out.println("avrHV-fun"+fun+"= "+sumHV/runtimes);
	 }//for-fun
  }//main
}
