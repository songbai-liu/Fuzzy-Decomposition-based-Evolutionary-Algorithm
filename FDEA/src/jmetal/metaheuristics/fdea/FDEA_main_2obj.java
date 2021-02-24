package jmetal.metaheuristics.fdea;

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
import jmetal.problems.MaF.MaF1;
import jmetal.problems.MaF.MaF10;
import jmetal.problems.MaF.MaF11;
import jmetal.problems.MaF.MaF12;
import jmetal.problems.MaF.MaF13;
import jmetal.problems.MaF.MaF2;
import jmetal.problems.MaF.MaF3;
import jmetal.problems.MaF.MaF4;
import jmetal.problems.MaF.MaF5_Concave;
import jmetal.problems.MaF.MaF6;
import jmetal.problems.MaF.MaF7;
import jmetal.problems.MaF.MaF8;
import jmetal.problems.MaF.MaF9;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG41;
import jmetal.problems.WFG.WFG42;
import jmetal.problems.WFG.WFG43;
import jmetal.problems.WFG.WFG44;
import jmetal.problems.WFG.WFG45;
import jmetal.problems.WFG.WFG46;
import jmetal.problems.WFG.WFG47;
import jmetal.problems.WFG.WFG48;
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

public class FDEA_main_2obj{
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
		
		for(int fun=13;fun<=16;fun++){
		int runtimes=1;
		double[] IGDarray=new double[runtimes];	
		double[] HVarray = new double[runtimes];
		long Execution_time=0;
		
		Problem problem = null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator\
		
		Operator mutation; // Mutation operator
		Operator selection; //Selection operator
		
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
			if(fun==6){
				//problem = new DTLZ1("Real",10,2);
			      problem = new DTLZ1("Real",7,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ1_3D.txt" );
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",12,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ2_3D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",12,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ2_3D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",12,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ2_3D.txt" );
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",12,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ5_3D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",12,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ5_3D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",22,3);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF\\DTLZ\\DTLZ7_3D.txt" );
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",2,20,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG\\WFG1_2D_5000.txt");
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",2,20,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG\\WFG2_2D_5000.txt");
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",2,20,2);
		  	      
		  	      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG\\WFG3_2D_5000.txt");
		  	    	}
			if(fun==16){
			      problem = new WFG9("Real",2,20,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG\\WFG9_2D_5000.txt");
			    	}//problem = new WFG1("Real");
			 if(fun==17){
		  	      problem = new WFG41("Real",2,20,2);
		  	      
		  	    indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG41_2D_5000.txt");
		  	    	}//problem = new WFG1("Real");
		  	if(fun==18){
		  	      problem = new WFG42("Real",2,20,2);
		  	      
		  	    indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG42_2D_5000.txt");
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
		  	      problem = new WFG43("Real",2,20,2);
		  	      
		  	    indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG43_2D_5000.txt");
		  	    	}
			if(fun==20){
			      problem = new WFG44("Real",2,20,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG44_2D_5000.txt");
			    	}//problem = new WFG1("Real");
			if(fun==21){
			      problem = new WFG45("Real",2,20,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG45_2D_5000.txt");
			    	}
		  	if(fun==22){
		  	      problem = new WFG46("Real",2,20,2);
		  	      
		  	    indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG46_2D_5000.txt");
		  	    	}//problem = new WFG1("Real");
		  	if(fun==23){
			      problem = new WFG47("Real",2,20,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG47_2D_5000.txt");
			    	}//problem = new WFG1("Real");
			if(fun==24){
			      problem = new WFG48("Real",2,20,2);
			      
			      indicators = new QualityIndicator(problem,"E:\\TruePF_Sampling\\PF\\WFG4X\\WFG48_2D_5000.txt");
			}
		} // else

		algorithm = new FDEA(problem);
		
		//algorithm = new NHaEA_Max(problem);
		
		algorithm.setInputParameter("maxGenerations", 300);
		algorithm.setInputParameter("populationSize", 100);
		algorithm.setInputParameter("T", 10);
		algorithm.setInputParameter("div1", 14);
		algorithm.setInputParameter("div2", 0);
		
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
		/*selection = SelectionFactory.getSelectionOperator("BinaryTournament",
				parameters);*/

		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);
		for(int i=0;i<runtimes;i++){
		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		Execution_time+=(System.currentTimeMillis() - initTime);
		//population.printObjectivesToFile("ADEA7_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+ (i+1)  + ".txt");
		/*if(fun>=29 && fun<=30){
			population.printVariablesToFile("NHaEA_Min_Variables_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"run_"+ (i+1)  + ".txt");
		}*/
		//IGDarray[i]=indicators.getIGD1(population);
		/*wfghvCalculator1 wfg = new wfghvCalculator1(population,fun);
		HVarray[i] = wfg.calculatewfghv();*/
		}
		//printGD("IGD"+ "-NSGAIII-"+problem.getName(),IGDarray);
		//printGD("NHaEA_Min_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_HV.txt",HVarray);
		double sumIGD=0;
		//double sumHV=0.0;
		for(int i=0;i<runtimes;i++){
			  sumIGD+=IGDarray[i];
			  //sumHV+=HVarray[i];
			}
		logger_.info("Total execution time: " + Execution_time + "ms");
	    System.out.println("avrIGD-fun"+fun+"= "+sumIGD/runtimes);
		//System.out.println("avrHV-fun"+fun+"= "+sumHV/runtimes);
	 }//for-fun
  }//main
}
