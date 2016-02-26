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

}	
