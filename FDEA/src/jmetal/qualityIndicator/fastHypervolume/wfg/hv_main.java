package jmetal.qualityIndicator.fastHypervolume.wfg;

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
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
import jmetal.problems.MaF.MaF1;
import jmetal.problems.MaF.MaF2;
import jmetal.problems.MaF.MaF3;
import jmetal.problems.MaF.MaF4;
import jmetal.problems.MaF.MaF5_Convex;
import jmetal.problems.MaF.MaF6;
import jmetal.problems.MaF.MaF7;
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
import jmetal.problems.cec2009Competition.UF10;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.util.Configuration;
import jmetal.util.JMException;

public class hv_main {
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
		
		for(int fun=13;fun<=21;fun++){
		int runtimes=30;
		double[] IGDarray=new double[runtimes];	
		double[] HVarray = new double[runtimes];
		double[] times = new double[runtimes];
		long Execution_time=0;
		
		Problem problem = null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
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
			if(fun==1){
		  	      problem = new ZDT1("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\ZDT1.pf" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==2){
		  	      problem = new ZDT2("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\ZDT2.pf" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==3){
		  	      problem = new ZDT3("Real");
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\ZDT3.pf" ) ;
		  	    	}
		  	if(fun==4){
			      problem = new ZDT4("Real",10);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\ZDT4.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==5){
			      problem = new ZDT6("Real",10);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\ZDT6.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==6){
				//problem = new DTLZ1("Real",10,2);
			      problem = new DTLZ1("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ1.pf") ;
			    	}
		  	if(fun==7){
			      //problem = new DTLZ2("Real",10,2);
		  		problem = new DTLZ2("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ2.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==8){
			      problem = new DTLZ3("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ3.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==9){
			      problem = new DTLZ4("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ4.pf" ) ;
			    	}
		  	if(fun==10){
			      problem = new DTLZ5("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ5.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==11){
			      problem = new DTLZ6("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ6.pf") ;
			    	}//problem = new WFG1("Real");
			if(fun==12){
			      problem = new DTLZ7("Real",10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\DTLZ7.pf" ) ;
			    	}
			if(fun==13){
		  	      problem = new WFG1("Real",4,10,3);
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG1.3D.pf" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==14){
		  	      problem = new WFG2("Real",4,10,3);
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG2.3D.pf" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==15){
		  	      problem = new WFG3("Real",4,10,3);
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG3.3D.pf" ) ;
		  	    	}
			if(fun==16){
			      problem = new WFG4("Real",4,10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG4.3D.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==17){
			      problem = new WFG5("Real",4,10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG5.3D.pf" ) ;
			    	}
		  	if(fun==18){
		  	      problem = new WFG6("Real",4,10,3);
		  	      
		  	      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG6.3D.pf" ) ;
		  	    	}//problem = new WFG1("Real");
		  	if(fun==19){
			      problem = new WFG7("Real",4,10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG7.3D.pf" ) ;
			    	}//problem = new WFG1("Real");
			if(fun==20){
			      problem = new WFG8("Real",4,10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG8.3D.pf" );
			    	}
			if(fun==21){
			      problem = new WFG9("Real",4,10,3);
			      
			      indicators = new QualityIndicator(problem,"H:\\truePF\\WFG9.3D.pf" ) ;
			    	}//problem = new WFG1("Real");
		  	if(fun==22){
			      problem = new MaF1("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			if(fun==23){
			      problem = new MaF2("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==24){
			      problem = new MaF3("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}//problem = new WFG1("Real");
			/*if(fun==25){
			      problem = new MaF4("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==26){
			      problem = new MaF5("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==27){
			      problem = new MaF6("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}
			if(fun==28){
			      problem = new MaF7("Real",12,3);
			      
			      indicators = new QualityIndicator(problem) ;
			    	}*/
			if(fun==25){
			      problem = new UF1("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF1.dat" ) ;
			    	}
			if(fun==26){
			      problem = new UF2("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF2.dat" ) ;
			    	}
			if(fun==27){
			      problem = new UF3("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF3.dat" ) ;
			    	}
			if(fun==28){
			      problem = new UF4("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF4.dat" ) ;
			    	}
			if(fun==29){
			      problem = new UF5("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF5.dat" ) ;
			    	}
			if(fun==30){
			      problem = new UF6("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF6.dat" ) ;
			    	}
			if(fun==31){
			      problem = new UF7("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF7.dat" ) ;
			    	}
			if(fun==32){
			      problem = new UF8("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF8.dat" ) ;
			    	}
			if(fun==33){
			      problem = new UF9("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF9.dat" ) ;
			    	}
			if(fun==34){
			      problem = new UF10("Real");
			      //indicators = new QualityIndicator(problem);
			      indicators = new QualityIndicator(problem,"H:\\truePF\\UF10.dat" ) ;
			    }
			if(fun==35){
				problem = new MOP1("Real",10);
			      
			   // indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP1.txt" ) ;
			}
          if(fun==36){
          	problem = new MOP2("Real",10);
			      
			   // indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP2.txt" ) ;
			}
          if(fun==37){
          	problem = new MOP3("Real",10);
			      
			    //indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP3.txt" ) ;
			}
          if(fun==38){
          	problem = new MOP4("Real",10);
			      
			    //indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP4.txt" ) ;
			}
          if(fun==39){
          	problem = new MOP5("Real",10);
			      
			    //indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP5.txt" ) ;
			}
          if(fun==40){
          	problem = new MOP6("Real",10,3);
			      
			   // indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP6.txt" ) ;
			}
          if(fun==41){
          	problem = new MOP7("Real",10,3);
			    
			    //indicators = new QualityIndicator(problem) ;
			    indicators = new QualityIndicator(problem,"H:\\truePF\\MOP7.txt" ) ;
			}
		} // else
		jmetal.qualityIndicator.util.MetricsUtil utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
		for(int i=0;i<runtimes;i++){
			String path = "E:\\MOEAFD\\MOEAD_ACD\\MOEAD_ACD_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_run"+(i+1)+".txt";
			double[][] population = utils_.readFront(path);
		    wfghvCalculator3 wfg = new wfghvCalculator3(population,fun,problem.getNumberOfObjectives());
		    HVarray[i] = wfg.calculatewfghv();
		}
		printGD("MOEAD_ACD_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+"_HV.txt",HVarray);
		double sumHV=0.0;
		for(int i=0;i<runtimes;i++){
			  sumHV+=HVarray[i];
			}
		logger_.info("Total execution time: " + Execution_time + "ms");
		System.out.println("avrHV-fun"+fun+"= "+sumHV/runtimes);
	 }//for-fun
  }//main
}
