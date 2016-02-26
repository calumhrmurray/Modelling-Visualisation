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
		// if current bonds negative total => positive deltaU & vice versa
   	        return 2.0 * J * l[i][j] * (lt + r + t + b);
   	} 
	

	static int[][] glauber(int[][] l, double temperature, final int W){
		// number of steps before displayed on screen
                for (int k=0; k<W*W; k++){ 
			// choose a random row and column
			int i = (int) (Math.random()*W); 
                	int j = (int) (Math.random()*W);
                	// work out energy change if spin is swapped
			double delta = deltaU(l,W,i,j);
			// metropolis	
			if (Math.random()<Math.exp(-delta/temperature*kb)){l[i][j] = -l[i][j];} 
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
			// make sure they're different, should I bother?
			while (i == s && j == t){
			      s = (int)(Math.random()*W); 
 			      t = (int)(Math.random()*W); 	
			}
			// Only perform calculation if the spins are not equal, stops unneccesary calc?
			if(l[i][j]!=l[s][t]){ 
				// work out energy change if spin is swapped 
				double deltaE = deltaU(l,W,i,j)+deltaU(l,W,s,t); 
				// fix for double counting
				if (i+1==s||i-1==s||j+1==t||j-1==t && deltaE>0){deltaE = deltaE - 2*J;} //should add +4*J with nearest neighbour, fix check for n.n. with periodic b.c.
				if (i+1==s||i-1==s||j+1==t||j-1==t && deltaE<=0){deltaE = deltaE + 2*J;} 
				// metropolis
				if (Math.random()<Math.exp(-deltaE/temperature*kb)){l[i][j]=-l[i][j];l[s][t]=-l[s][t];}
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
		// Chooses simulation method
		if(sim==0){
		Animate.main(args);
		} else if (sim==1) { Simulate.startMeasurements(l,d,W,1000,J,kb,output); 
		} else {System.out.println("Simulation does not exist: sim = 0 for animation, sim = 1 for Chi and Cv calculations");}
			

		}

}

	//----------------------------------------------------------------------


