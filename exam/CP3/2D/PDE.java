/*
 *	PDE solver
 */
import java.io.*;

public class PDE{

	private static double[][] c;
	private static double[][] cPrime;
	private static double[][] p;
	private static double deltaX = 1;
	private static double perm = 1;

	double[][] getPhiArray() {
		return c;
	}

	double[][] getPhiPrimeArray() {
		return cPrime;
	}

	//----------------------------------------------------------------------

	static double[][] initChargeDensity(final int W){
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ if(i==50 && j ==50){p[i][j] = 255;} else {p[i][j] = 0;}}}
		return p;
	}

	static double[][] init(final int W) {
		c = p;
		return c;
	}

	static double[][] update(final int W){
		jacobi(c,cPrime,W);
		c = cPrime;
		return c;
	}

	//----------------------------------------------------------------------

	static double[][] jacobi(double[][] array, double[][] arrayPrime, final int W){
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ if(i==W-1 && j==W-1 ){arrayPrime[i][j] = 0.25*(array[i-1][j]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(i==W-1 && j==0){arrayPrime[i][j] = 0.25*(array[i-1][j]+array[i][j+1]+deltaX*deltaX*p[i][j]);
					}else if(i==0 && j==W-1){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(i==0 && j==0){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i][j+1]+deltaX*deltaX*p[i][j]);
					}else if(i==W-1){arrayPrime[i][j] = 0.25*(array[i-1][j]+array[i][j+1]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(j==W-1){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i-1][j]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(j==0){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i-1][j]+array[i][j+1]+deltaX*deltaX*p[i][j]);
					}else if(i==0){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i][j+1]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else{arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i-1][j]+array[i][j+1]+array[i][j-1]+deltaX*deltaX*p[i][j]);}}}
		return arrayPrime;
	}



	//----------------------------------------------------------------------


	public static void main(final String[] args) throws Exception {
		if (args.length != 3) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final int acc = Integer.parseInt(args[1]);
		final int out = Integer.parseInt(args[2]);

		c = new double[W][W];
		cPrime = new double[W][W];
		p = new double[W][W];		

		// Initialise model
		initChargeDensity(W);
		init(W);

		Animate.main(args);
		//Printer.main(args); 

	}

}

	//----------------------------------------------------------------------

/**


		What's supposed to happen at the boundary all zero at the edge or beyond?
		Think of a better boundary sln (might be in your notes)
	

















**/	
