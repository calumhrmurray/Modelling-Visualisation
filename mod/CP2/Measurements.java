/**
 *  Contains different simulations to be run by Ising , simplifies input 
 *
 */
import java.io.*;

public class Measurements{

	public static void phaseDiagram(int[][] l,final int W, final double p2 ,PrintWriter output){

		// initialise doubles and ints for measurements
		double I; double avFracI; double avI; double varI;
		int ndata;
		// number of measurements taken
		int N = 30;
		// number of steps waiting for equilibrium
		int wait = 10; 
		double[] sArray = new double[N];

		// for measurements
		// loop over a range of p1,p3
		for (double p1=0;p1<1;p1+=0.01){
		for (double p3=0;p3<1;p3+=0.01){
			// needs to be initialised for every probability
			l = Sirs.init(W);
                        // set measurements to zero each time around
			I = 0; avI =0; varI = 0;			
			ndata = 0;
	
			//perform measurements
			for (int i=0;i<(wait+1+N*10);i++){
				l = Sirs.update(W,p1,p2,p3);		
				// wait till equilibrium is reached
				if(i>wait){
					// take measurements every 10 steps
					if(i % 10 == 0){	
						// take measurements and store them
						sArray[ndata] = updateI(l,W);
						ndata++;						
					}
				}
			}
			// calculate Cv,Chi and their errors
			avI = getMean(sArray,ndata);
			varI = getVariance(sArray,ndata);	
			// output stuff
	                output.printf(" "+p1+" "+p3+" "+avI+" "+varI+"\n");
			      if (0.0999 < p1 && p1 < 0.1001 && p3 == 0
			        ||0.1999 < p1 && p1 < 0.2001 && p3 ==0
				||0.2999 < p1 && p1 < 0.3001 && p3 ==0
				||0.3999 < p1 && p1 < 0.4001 && p3 ==0
				||0.4999 < p1 && p1 < 0.5001 && p3 ==0
				||0.5999 < p1 && p1 < 0.6001 && p3 ==0
				||0.6999 < p1 && p1 < 0.7001 && p3 ==0
				||0.7999 < p1 && p1 < 0.8001 && p3 ==0
				||0.8999 < p1 && p1 < 0.9001 && p3 ==0){
			System.out.println("Progress:" +p1); // nothing in array
			}
			if(0.9999 <p1 && p1 < 1.0001 && 0.9999 <p3 && p3 < 1.0001){System.out.println("Complete");}

		}
		// new line for new x values
		output.printf("\n");
		}		
		output.close();	

	} 

	public static void infectedTimeDiagram(int[][] l, final int W, final double p1, final double p2, final double p3,PrintWriter output){
	double I;
	// plot I/N over time
	for (int t=0; t<10001; t++){
			I = 0;
			// calculate number of infected lattice sites
			I = updateI(l,W);
			// update the lattice (perform 1 MC sweep)
			l = Sirs.update(W,p1,p2,p3);		
			// output the number of infected sites per lattice site against time
		        output.printf(" "+t+ " "+I+"\n");
			if ( t % 100 == 0){System.out.printf("\r Time: "+t);}
		}
		output.close();	

	}

	// for each (p1,im) measure the fraction infected after some time
	// graph infected(p1)
	public static void immunityDiagram(int[][] l, final int W, final double p2, final double p3, PrintWriter output){

		// number of measurements taken for each data point
		int N = 30;
		// number of steps waiting for equilibrium
		int wait = 10; 
	        // array to store measurements
		double[] sArray = new double[N];

		for (double im = 0; im<1; im+=0.1){
		for (double p1 = 0; p1<1; p1+=0.1){

			l = Sirs.initImmune(W,im);

	        	// initialise doubles and ints for measurements
			double I; double avI; double varI;
			int ndata =0;	
		        
			//perform measurements
			for (int i=0;i<(wait+1+N*10);i++){
				l = Sirs.update(W,p1,p2,p3);		
				// wait till equilibrium is reached
				if(i>wait){
					// take measurements every 10 steps
					if(i % 10 == 0){	
						// take measurements and store them
						sArray[ndata] = updateI(l,W);
						ndata++;						
					}
				}
			}
			// calculate Cv,Chi and their errors
			avI = getMean(sArray,ndata);
			varI = getVariance(sArray,ndata);	
			// output stuff
	                output.printf(" "+im+" "+p1+" "+avI+" "+varI+"\n");
			      if (0.0999 < p1 && p1 < 0.1001 && im == 0
			        ||0.1999 < p1 && p1 < 0.2001 && im ==0
				||0.2999 < p1 && p1 < 0.3001 && im ==0
				||0.3999 < p1 && p1 < 0.4001 && im ==0
				||0.4999 < p1 && p1 < 0.5001 && im ==0
				||0.5999 < p1 && p1 < 0.6001 && im ==0
				||0.6999 < p1 && p1 < 0.7001 && im ==0
				||0.7999 < p1 && p1 < 0.8001 && im ==0
				||0.8999 < p1 && p1 < 0.9001 && im ==0){
			System.out.println("\r Progress:" +p1); // nothing in array
			}
			if(0.9999 <p1 && p1 < 1.0001 && im == 0){System.out.println("\r Almost there...");}	
			}
		// new line for new x values
		output.printf("\n");
		}	

		// loop over a range of p1
		
		output.close();	

	}

	//----------------------------------------------------------------------
	// methods acting on lattice
	

	// returns the number of infected sites
	public static double updateI(int[][] l, final int W){
		double I = 0;
		// loop over all lattice points
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
			//
			if( l[i][j] == -1){ I++;}
		}}	
		return I/(W*W);
	}

	// returns the number of recovered sites
	public static double updateR(int[][] l, final int W){
		double I = 0;
		// loop over all lattice points
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
			//
			if( l[i][j] == 0){ I++;}
		}}	
		return I/(W*W);
	}

	// returns the number of sus sites
	public static double updateS(int[][] l, final int W){
		double I = 0;
		// loop over all lattice points
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
			//
			if( l[i][j] == 1){ I++;}
		}}	
		return I/(W*W);
	}
 
	//----------------------------------------------------------------------
	// stats methods


	//calculates average of an array of measurements
	public static double getMean(double[] data, int N){
		double average = 0;
		for (int i=0; i < N; i++){
			average += data[i];			
		}
		average = average/N;
		return average;
	}

		// Calculates variance of two quantities
	public static double getVariance(double[] data, int N){
		double variance = 0;
		double t0 = 0; double t02 = 0; 
		// sum up all the elements of each row
		for (int i = 0; i<N; i++){
			t0 += data[i];
			t02 += data[i]*data[i];
		}
		variance = t02/N-(t0/N)*(t0/N);
		return variance;
	}

	//----------------------------------------------------------------------

}