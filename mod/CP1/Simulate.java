/**
 *  Contains different simulations to be run by Ising , simplifies input 
 *
 */
import java.io.*;

public class Simulate{

	public static void startMeasurements(int[][] l, final int d, final int W, int N, double J, double kb,PrintWriter output){

		// initialise doubles and ints for measurements
		double M,Temp,E;
		int ndata;
		// number of steps waiting for equilibrium
		int wait = 100; 
		double[][] sArray = new double[N][2];
		// data array
		double[][] data = new double[1][4];

		// for measurements
		// loop over a range of T
		for (double T=12;T<50;T++){
                        // set measurements to zero each time around
			M = 0;
			E = 0;
			ndata = 0;
			Temp = T/10;		
			//perform measurements
			for (int i=0;i<(wait+1+N*10);i++){
				l = Ising.update(d,Temp,W);		
				// wait till equilibrium is reached
				if(i>wait){
					// take measurements every 10 steps
					if(i % 10 == 0){	
						// take measurements and store them
						sArray[ndata][0] = Stats.updateM(l,W);
						sArray[ndata][1] = Stats.updateE(l,W,J);
						ndata++;						
					}
				}
			}
			// calculate Cv,Chi and their errors
			data = Stats.stats(sArray, N, Temp, W, kb);
			//M = Stats.updateM(l,W);
			//E = Stats.updateE(l,W,J);
			// output stuff, ignore first measurement
			if ( T != 12){output.printf(" "+Temp+" "+data[0][0]+" "+data[0][1]+ " "+data[0][2]+" "+data[0][3]+" "+M+" "+E+" "+ndata+"\n");}
		}		
		output.close();	

	}

	public static void pottsMeasurements(int[][] l, final int d, final int W, int c, double J, double kb,PrintWriter output){

		int wait = 3000;
		double avE = 0;
		double avE2 = 0;
		double variance = 0;
		double error1 = 0;
		double error2 = 0;

		for (double temp = 0.3; temp < 2.6; temp += 0.1){

			l = Ising.initPotts(W,c);

			for ( int t = 0; t < wait+1000; t++){

				if (t>wait){
					avE += Stats.pottsUpdateE(l,W,J);
					avE2 = avE*avE;
				}

			l = Ising.updatePotts(d,temp,W);		
			System.out.printf("\r Progress:" +temp);
			}
	
			avE = avE/1000;
			avE2 = avE2/1000;
			variance = avE2-avE*avE;
			output.printf(" "+temp+" "+avE+" "+variance+ " "+error1+" "+error2+"\n");
			avE = avE2 = variance = 0;

		}


		output.close();	

	}

}	
