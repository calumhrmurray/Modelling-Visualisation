/**
 *	Static methods used to make measurements
 */

public class Stats{

	// returns the magnetisation of the lattice
	public static double updateM(int[][] l, final int W){
		double M = 0;
		// Sum all the spins
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
		M += l[i][j];
		}}	
		return M;
	}
	// returns the energy of the lattice
	public static double updateE(int[][] l, final int W, double J){
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
	// returns Chi Cv and errors in measurements
	public static double[][] stats(double[][] sArray, int N, double T, final int W, double kb){
		double[][] rArray = new double[N][2];
		double[][] dataArray = new double[1][2];
		double[][] sigArray = new double[1][2];
		double[][] statArray = new double[1][4];
		double tM = 0; double tE = 0;	double tM2 = 0; double tE2 = 0;
		// Chi and Cv
		dataArray = calc(sArray,N,T,W,kb);
		// Errors
		sigArray = bootstrap(sArray,N,T,W,kb);
		//sigArray = jacknife(sArray,dataArray,N,T,W,kb);
		// Chi, Cv				
		statArray[0][0] = dataArray[0][0];
		statArray[0][1] = dataArray[0][1];
		// errors
		statArray[0][2] = sigArray[0][0];
		statArray[0][3] = sigArray[0][1];
		return statArray;
	}
	// Only used within stats
	//----------------------------------------------------------------------
	// Calculator finds Cv and Chi
	public static double[][] calc(double[][] sArray, int N, double T, final int W, double kb){
		double[][] dataArray = new double[1][2];
		dataArray = variance(sArray,N,N);
		// Chi, Cv				
		dataArray[0][0] = (dataArray[0][0])/(kb*T*W*W);
		dataArray[0][1] = (dataArray[0][1])/(kb*T*T);
		return dataArray;
	}
	// bootstrap
	public static double[][] bootstrap(double[][] sArray, int N, double T, final int W, double kb){
		// number of Cv,Chi taken from bootstrap
		int n = 10;
		double[][] rArray = new double[N][2];
		double[][] bootArray = new double[n][2];
		double[][] sigArray = new double[1][2];
		// calc. Cv and Chi n times using resampled data
		for (int i = 0; i<n; i++){
			// pick out n measurements at random and compute Cv, Chi
			for (int j = 0; j<N; j++){
				rArray[j][0] = sArray[(int)(Math.random()*N)][0];	
				rArray[j][1] = sArray[(int)(Math.random()*N)][1];	
			}
			// array with n values of Cv and Chi
			bootArray[i][0] = calc(rArray,N,T,W,kb)[0][0];
			bootArray[i][1] = calc(rArray,N,T,W,kb)[0][1];
		}	
		// calculate sum of (ci-c)2
		sigArray = variance(bootArray,n,N);
		sigArray[0][0] = Math.sqrt(sigArray[0][0]);
		sigArray[0][1] = Math.sqrt(sigArray[0][1]);
		return sigArray	;
	}
	// jacknife
	public static double[][] jacknife(double[][] sArray,double[][] dataArray, int N, double T, final int W, double kb){
		// number of Cv,Chi taken from bootstrap
		int r;
		int n = 10;
		double chiError = 0; double CvError = 0;
		double[][] jArray = new double[N-1][2];
		double[][] sigArray = new double[1][2];
		// calc. Cv and Chi n times using resampled data
		for (int i = 0; i<n; i++){
			// remove the jth measurement and compute Cv, Chi
			for (int j = 0; j<N-1; j++){
   	                       if (j == N-1){ r = 0;} else{ r = j+1;}
				jArray[j][0] = sArray[r][0];	
				jArray[j][1] = sArray[r][1];	
			}
			// array with n values of Cv and Chi
			chiError += (calc(jArray,N-1,T,W,kb)[0][0]-dataArray[0][0])*(calc(jArray,N-1,T,W,kb)[0][0]-dataArray[0][0]);
			CvError += (calc(jArray,N-1,T,W,kb)[0][1]-dataArray[0][1])*(calc(jArray,N-1,T,W,kb)[0][1]-dataArray[0][1]);			
		}			
		sigArray[0][0] = Math.sqrt(chiError);
		sigArray[0][1] = Math.sqrt(CvError);
		return sigArray	;
	}
	// Calculates variance of two quantities
	public static double[][] variance(double[][] data, int n, int N){
		double[][] output = new double[1][2];
		double t0 = 0; double t1 = 0; double t02 = 0; double t12 = 0;
		// sum up all the elements of each row
		for (int i = 0; i<n; i++){
			t0 += data[i][0];
			t1 += data[i][1];
			t02 += data[i][0]*data[i][0];
			t12 += data[i][1]*data[i][1];
		}

		output[0][0] = t02/N-(t0/N)*(t0/N);
		output[0][1] = t12/N-(t1/N)*(t1/N);

		return output;
	}


}
