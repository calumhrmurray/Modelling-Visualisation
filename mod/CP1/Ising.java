/**
 * Methods etc. used in the Ising model, also runs animate class
 */
import java.lang.Math;
import java.io.*;

public class Ising{

	// Constants
	static final double J = 1;
	static final double kb = 1;

	private static int[][] l;

	int[][] getArray() {
		return l;
	}

	//----------------------------------------------------------------------

	static int[][] update(final int d, double temp, int W){
			if(d == 0){	
			l = glauber(l,temp,W);
			} else {
			l = kawasaki(l,temp,W);
			}	
         return l;
	}
	static int[][] init(final int W, final int c) {
		if (c == 0){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){l[i][j] = 1;}}
        	} else if (c==1){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){l[i][j] = -1;}}
		} else if (c==2){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){if (Math.random() < 0.5){ l[i][j] = 1;} else{ l[i][j] = -1;}}}
		} else {
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){if (i < W/2){ l[i][j] = 1;} else{ l[i][j] = -1;}}}
		}	
         return l;	
	}	

	//----------------------------------------------------------------------


	// Compute energy change from changeing the spin of a lattice site
	static double deltaU(int[][] l,final int W,int i, int j) {
		int lt,r,t,b;  // values of neighboring spins
		// act on the array l
   	        if (i == 0) lt = l[W-1][j]; else lt = l[i-1][j];
   	        if (i == W-1) r = l[0][j]; else r = l[i+1][j];
   	        if (j == 0) t = l[i][W-1]; else t = l[i][j-1];
   	        if (j == W-1) b = l[i][0]; else b = l[i][j+1];
   	        return 2.0 * l[i][j] * (lt + r + t + b);
   	} 
	

	static int[][] glauber(int[][] l, double temperature, final int W){
		// number of steps before displayed on screen
                for (int k=0; k<W*W; k++){ 
			// choose a random row and column
			int i = (int) (Math.random()*W); 
                	int j = (int) (Math.random()*W);
                	// work out energy change if spin is swapped
			double delta = J*deltaU(l,W,i,j);
			// metropolis	
			if (Math.random()<Math.min(1,Math.exp(-delta/temperature*kb))){l[i][j] = -l[i][j];} 
		}			
		return l;
		}	


	static int[][] kawasaki(int[][] l, double temperature, final int W){
 		// number of steps before displayed on screen
		for (int k=0; k<W*W; k++){ 
			// choose a random row and column
			int i = (int)(Math.random()*W); int j = (int)(Math.random()*W);
			// choose a new random row and column  
			int s = (int)(Math.random()*W); int t = (int)(Math.random()*W); 
			// make sure they're different, should I bother????????????
			while (i == s && j == t){
			      s = (int)(Math.random()*W); 
 			      t = (int)(Math.random()*W); 	
			}
			// Only perform calculation if the spins are not equal, stops unneccesary calc?
			if(l[i][j]!=l[s][t]){ 
				// work out energy change if spin is swapped 
				double deltaE = deltaU(l,W,i,j)+deltaU(l,W,s,t); 
				// fix for double counting
				if (i+1==s||i-1==s||j+1==t||j-1==t && deltaE>0){deltaE = deltaE - 2;} 
				// metropolis
				// Will math.random return numbers only less than 1?		
				// Why do min of this????????? surely works anyway???????????	
				if (Math.random()<Math.min(1,Math.exp(-deltaE/temperature*kb))){l[i][j] = -l[i][j];l[s][t] = -l[s][t];}
			}
		}			
		return l;
		}
	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 5) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final int t = Integer.parseInt	(args[1]);
		final int c = Integer.parseInt(args[2]);
		final int d = Integer.parseInt(args[3]);
		final int sim = Integer.parseInt(args[4]);

		String filename = "random.data";
		// Set output filename
		if(d==0){filename = "glauber.data";}else{filename = "kawasaki.data";}
		// start printer
		PrintWriter output = new PrintWriter(new FileWriter(filename));

		l = new int[W][W];
		// Initialise model	
		init(W,c);

		if(sim==0){
		Animate.main(args);
		} else if (sim==1) { Simulate.startMeasurements(l,d,W,1000,J,kb,output); 
		} else {System.out.println("Simulation does not exist: sim = 0 for animation, sim = 1 for Chi and Cv ");}
			

		}





}

	//----------------------------------------------------------------------


	//** Check list **//

	/**	
					

		\/\/\/\/\/Sort out kawasaki - definetly incorrect!!!!\/\/\/\/\/
		\/\/\/\/\/Neaten up code\/\/\/\/\/\/
   	        \/\/\/\/\/Include J and k \/\/\/\/\/
		\/\/\/\/\/nicer measurements methods\/\/\/\/
		\/\/\/\/\/make Cv measurements\/\/\/\/\/
		\/\/\/\/\/work out how to do errors\/\/\/\/
		\/\/\/\/\/sample every 10 steps\/\/\/\/\/
		\/\/\/\/\/fix chi\/\/\/\/\/
		\/\/\/\/\/neaten up errors\/\/\/\/
		\/\/\/\/\/combine errors with Chi and Cv?\/\/\/\/
		go through checkpoint get all graphs ready
		\/\/\/\/\/GENERALISE methods\/\/\/\/\/
		\/\/\/\/\/\/could add an rms method\/\/\/\/
		\/\/\/\/\/\/\could bootstrap method be simplified?\/\/\/\/\/
		\/\/\/\/\//\both glauber and kawasaki animate seems to hit critical T at the wrong value! SEEMS TO WORK NOW\/\/\/\/
		check that kawasaki double counting is defo correct
		\/\/\/\/\/\/\/add jacknife\/\/\/\/\/\/
		sort out the wrong first value
		add text block with total magnets and energy to animate output
		add start button
		fix low T kawasaki



	*/


/*		// number of steps waiting for equilibrium
		int wait = 100; 
		// number of measurements taken for each value of T
		int N = 1000;	
		// initialise doubles and ints for measurements
		double M,Temp,E;
		int ndata;
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
				l = update(d,Temp,W);		
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
			M = Stats.updateM(l,W);
			E = Stats.updateE(l,W,J);
			// output stuff
	                output.printf(" "+Temp+" "+data[0][0]+" "+data[0][1]+ " "+data[0][2]+" "+data[0][3]+" "+M+" "+E+" "+ndata+"\n");
		}		
		output.close();	
			}	
*/



