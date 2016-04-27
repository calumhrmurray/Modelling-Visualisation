/**
 * Methods etc. used in the Ising model, also runs animate class
 */
import java.lang.Math;
import java.io.*;

public class Super{

	// Constants
	static final double J = 1;
	static final double kb = 1;

	private static int[][] l;

	int[][] getArray() {
		return l;
	}

	//----------------------------------------------------------------------

	static int[][] update(double T,int W){
		if(Math.random()<0.5){
			l = kawasaki(l,T,W);
		} else {
			l = glauber(l,T,W);
		}
         return l;
	}

	static int[][] init(final int W, final double c) {
		for (int i=0; i<W; i++){ 
		for (int j=0; j<W; j++){ 
		if (Math.random() < c){
			l[i][j] = 0;
		} else {
			if (Math.random() < 0.5){l[i][j]=-1;}else{l[i][j]=1;}
		} 
		}}	
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
			// make sure they're different and not next to one another
			while (i == s && j == t || (i+1==W?0:i+1)==s ||(i-1==0?W-1:i-1)==s || (j+1==W?0:j+1)==t || (j-1==0?W-1:j-1)==t ){
			      s = (int)(Math.random()*W);  
 			      t = (int)(Math.random()*W); 	
			}
			// Only perform calculation if the spins are not equal, stops unneccesary calc?
			if(l[i][j]!=l[s][t]){ 
				// work out energy change if spin is swapped 
				double deltaE = deltaU(l,W,i,j)+deltaU(l,W,s,t); 
				// metropolis
				int spot = l[i][j];
				if (Math.random()<Math.exp(-deltaE/temperature*kb)){l[i][j]=l[s][t];l[s][t]=spot;}
			} 
		}			
		return l;
		}

		
	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 4) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]); // lattice size
		final double T = Double.parseDouble(args[1]); // initial temperature
		final double c = Double.parseDouble(args[2]); // initialisation
		final int sim = Integer.parseInt(args[3]); // simulation

		String filename = "random.data";
		// Set output filename
		filename = "SuperFluid.data";
		// start printer
		PrintWriter output = new PrintWriter(new FileWriter(filename));

		l = new int[W][W];
		// Initialise model	
		init(W,c);

		if (sim==0){SuperAnimate.main(args);}else{startMeasurements(l,W,output);}

		}

	//----------------------------------------------------------------------

	// returns the energy of the lattice
	public static double updateE(int[][] l, final int W){
		int r,t;
		double E =0;
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
   	        if (i == W-1){r = l[0][j];}else{r = l[i+1][j];}
   	        if (j == 0){t = l[i][W-1];}else{t = l[i][j-1];}
		E += -J*l[i][j]*(r+t);
		}}
		return E;
	}

	// single bootstrap N is the number of measurements in the array
	public static double bootstrap(double[] array, int N, double T){
		// number of Cv calculated
		int n = 10;
		double eStore = 0; double E = 0; double E2 = 0;
		double[] recalculated = new double[n];
		double error = 0;
		double avC = 0;
		double avC2 = 0;
		// calc. c n times using resampled data
		for (int i = 0; i<n; i++){
			// pick out n measurements at random and compute Cv, Chi
			for (int j = 0; j<N; j++){
				eStore = array[(int)(Math.random()*N)];
				E += eStore;
				E2 += eStore*eStore;	
			}
			// array with n values of Cv 
			E = E/N;
			E2 = E2/N;			
			recalculated[i] = (E2-E*E)/(kb*T*T);
		}	
		// taking the recalculated array we calculate sigma on this recalculated array
		for (int l = 0; l<n; l++){
			avC += recalculated[l];
			avC2 += recalculated[l]*recalculated[l];
		}
		avC = avC/n;
		avC2 = avC2/n;
		error = Math.sqrt(avC2-avC*avC);
		return error;
	}

	public static void startMeasurements(int[][] l, final int W,PrintWriter output){

		double E = 0;
		double E2 = 0;
		double Cv = 0;
		int N = 500;
		int wait = 1000;
		double[] eStore = new double[N];	
		int ndata = 0;
		double error = 0;
		int step = 5;

		// for measurements
		// loop over a range of T
		for (double T=0.1;T<1.51;T+=0.04){
			//reinit lattice
			init(W,0.8);
                        // set measurements to zero each time around
			E = 0;
			E2 = 0;
			ndata = 0;
			//perform measurements
			for (int i=0;i<(wait+1+N*step);i++){
				l = Super.update(T,W);		
				// wait till equilibrium is reached
				if(i>wait){
					// take measurements every step steps
					if(i % step == 0){	
						// take measurements and store them
						eStore[ndata] = updateE(l,W);
						E += eStore[ndata];
						E2 += eStore[ndata]*eStore[ndata];
						ndata++;
					}
				}
			}
			// calculate Cv
			E = E/N;
			E2 = E2/N;
			Cv = (E2-E*E)/(kb*T*T);
			error = bootstrap(eStore,N,T);
			output.printf(" "+T+" "+E+" "+Cv+" "+error+"\n");
			System.out.printf("\r Progress:" +T);
		}	
		System.out.printf("\n");	
		output.close();	

	}

}

