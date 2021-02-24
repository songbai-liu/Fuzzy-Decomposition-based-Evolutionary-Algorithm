//  MOEAD_main.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.moead;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
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
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
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
import jmetal.problems.cec2009Competition.*;
import jmetal.problems.mDTLZ.mDTLZ1;
import jmetal.problems.mDTLZ.mDTLZ2;
import jmetal.problems.mDTLZ.mDTLZ3;
import jmetal.problems.mDTLZ.mDTLZ4;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator3;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

/**
 * This class executes the algorithm described in: H. Li and Q. Zhang,
 * "Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D
 * and NSGA-II". IEEE Trans on Evolutionary Computation, vol. 12, no 2, pp
 * 284-302, April/2009.
 */
public class computeIGD_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	
	public static jmetal.qualityIndicator.util.MetricsUtil utils_;

	/**
	 * @param args
	 *            Command line arguments. The first (optional) argument
	 *            specifies the problem to solve.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options - jmetal.metaheuristics.moead.MOEAD_main
	 *             - jmetal.metaheuristics.moead.MOEAD_main problemName -
	 *             jmetal.metaheuristics.moead.MOEAD_main problemName
	 *             ParetoFrontFile
	 * @throws ClassNotFoundException
	 */
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	      BufferedWriter bw      = new BufferedWriter(osw)        ;               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");
	        bw.newLine();        
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		
		// Logger object and file to store log messages
				logger_ = Configuration.logger_;
				fileHandler_ = new FileHandler("MOEAD.log");
				logger_.addHandler(fileHandler_);
		for(int fun=6;fun<=26;fun++){	
		int runtimes=30;
	
		double[] IGDarray=new double[runtimes];
		double[] HVarray=new double[runtimes];
		long Execuion_time=0;
		
		Problem problem=null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		QualityIndicator indicators; // Object to get quality indicators

		HashMap parameters; // Operator parameters

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
			      problem = new DTLZ1("Real",14,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ1_10D.txt" );
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ2_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ2_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ2_10D.txt" );
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ5_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",19,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ5_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",29,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ7_10D.txt" );
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",18,20,10);
		  	      
		  	      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG1_10D.txt" );
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",18,20,10);
		  	      
		  	      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG2_10D.txt" );
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",18,20,10);
		  	      
		  	      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG3_10D.txt" );
		  	    	}
			if(fun==16){
			      problem = new WFG4("Real",18,20,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG4_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==17){
			      problem = new WFG5("Real",18,20,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG4_10D.txt" );
			    	}
		  	if(fun==18){
		  	      problem = new WFG6("Real",18,20,10);
		  	      
		  	      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG4_10D.txt" );
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
			      problem = new WFG7("Real",18,20,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG4_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==20){
			      problem = new WFG8("Real",18,20,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG4_10D.txt" );
			    	}
			if(fun==21){
			      problem = new WFG9("Real",18,20,10);
			      
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\WFG\\WFG4_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==22){
			      problem = new MaF1("Real",19,10);
			      
			      //indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\MaF\\MaF1_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==23){
			      problem = new MaF2("Real",19,10);
			      
			     // indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\MaF\\MaF2_10D.txt" );
			    	}
			if(fun==24){
			      problem = new MaF3("Real",19,10);
			      
			      //indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\MaF\\MaF3_10D.txt" );
			    	}//problem = new WFG1("Real");
			if(fun==25){
			      problem = new MaF4("Real",19,10);
			      
			      //indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\MaF\\MaF4_10D.txt" );
			    	}
			if(fun==26){
			      //problem = new MaF5_Convex("Real",14,5);
			      problem = new MaF5_Concave("Real",19,10);
			      
			      //indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\MaF\\MaF5_10D.txt" );
			    	}
			if(fun==27){
			      problem = new MaF6("Real",19,10);
			      
			      //indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ5_10D.txt" );
			    	}
			if(fun==28){
			      problem = new MaF7("Real",29,10);
			      
			      //indicators = new QualityIndicator(problem) ;
			      indicators = new QualityIndicator(problem,"D:\\TruePF\\DTLZ\\DTLZ7_10D.txt" );
			    	}
			if(fun==29){
			      problem = new MaF8("Real",2,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==30){
			      problem = new MaF9("Real",2,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==31){
			      problem = new MaF10("Real",8,20,5);//WFG1
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==32){
			      problem = new MaF11("Real",8,20,5);//WFG2
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==33){
			      problem = new MaF12("Real",8,20,5);//WFG9
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==34){
			      problem = new MaF13("Real",5,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			 if(fun==35){
			      problem = new mDTLZ1("Real",15,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==36){
			      problem = new mDTLZ2("Real",15,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==37){
			      problem = new mDTLZ3("Real",15,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
	          if(fun==38){
			      problem = new mDTLZ4("Real",15,5);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
		} // else
		utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	
		
		for(int i=0;i<runtimes;i++){
	       //String path = "D:\\Research Achievement\\MOEA_AC_2019\\Running_data\\NSGA_III"+"\\NSGAIII_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
	       //String path = "D:\\Research Achievement\\MOEA_AC_2019\\Running_data\\TDEA"+"\\TDEA_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
	       //String path = "D:\\Research Achievement\\MOEA_AC_2019\\Running_data\\VaEA"+"\\VaEA_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
	       String path = "D:\\Research Achievement\\MOEA_AC_2019\\Running_data\\MaOEA_C"+"\\MaOEAC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
	       //String path = "D:\\Research Achievement\\MOEA_AC_2019\\Running_data\\MOEA_AD1"+"\\MOEAAD1_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
	       //String path = "D:\\Research Achievement\\MOEA_AC_2019\\Running_data\\MOEA_AD2"+"\\MOEAAD2_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
	       double[][] population = utils_.readFront(path);
	       IGDarray[i] = indicators.getIGD1(population);
	       wfghvCalculator3 wfg = new wfghvCalculator3(population,fun, problem.getNumberOfObjectives());
		   HVarray[i] = wfg.calculatewfghv();
		}
		printGD("MaOEAC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_IGD.txt",IGDarray);
		printGD("MaOEAC_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_HV.txt",HVarray);
	
		double sumIGD=0;
		double sumHV = 0;
		for(int i=0;i<runtimes;i++){
		  sumIGD+=IGDarray[i];
		  sumHV+=HVarray[i];
		}
		logger_.info("Total execution time: " + Execuion_time + "ms");
		System.out.println("avrIGD-fun"+fun+"= "+sumIGD/runtimes);
		System.out.println("avrHV-fun"+fun+" = "+sumHV/runtimes);
	 }//for
	} // main
} // MOEAD_main
