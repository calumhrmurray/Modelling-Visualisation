/*
 *	PDE solver
 */
import java.io.*;

public class PDE{

	private static double[][][] c;
	private static double[][][] p;
	private static double deltaX = 1;
	private static double perm = 1;

	double[][][] getPhiArray() {
		return c;
	}

	//----------------------------------------------------------------------

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] pointCharge(final int W){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}}
		p[(int)(W/2 + 0.5d)][(int)(W/2 + 0.5d)][(int)(W/2 + 0.5d)] = 100;
		return p;
	}

	static double[][][] init() {
		c = p;
		return c;
	}

	static double[][][] update(final int W){
		//for(int i = 0; i<1000; i++){
		return jacobi(W);
		//}
	}

	//----------------------------------------------------------------------

	static double[][][] jacobi(final int W){
		double[][][] arrayPrime = c;
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
		for (int k=0; k<W; k++){arrayPrime[i][j][k] = (c[i+1][j][k]+c[i==0?W:i-1][j][k]   // if an array element on the boundary is called value will be zero
									+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;}}}
		return c = arrayPrime;
	}
	//l[x==W-1?0:(x+1)][y]==-1

	// takes in the array filled with potential outputs magnitude of the eField and the direction
	static double[] eField(double[][][] array, int i, int j, int k, int W){
      		double[] field = new double[2];
		double x = 0; double y = 0;
		// with boundary conditions
		x = -(array[i+1][j][k]+array[i==0?W:i-1][j][k]-2*array[i][j][k])/(deltaX*deltaX);
		y = -(array[i][j+1][k]+array[i][j==0?W:j-1][k]-2*array[i][j][k])/(deltaX*deltaX);
		// magnitude of the field
		// all need to be between 0 and 1
		field[0] = Math.sqrt(x*x+y*y);
		// check this is in radians
		field[1] = Math.asin(y/(x*x+y*y));
		return field;			
		
	}

	//----------------------------------------------------------------------


	public static void main(final String[] args) throws Exception {
		if (args.length != 3) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final int acc = Integer.parseInt(args[1]);
		final int out = Integer.parseInt(args[2]);

		c = new double[W+1][W+1][W+1];
		p = new double[W+1][W+1][W+1];		

		// Initialise model
		p = pointCharge(W);
		c = p;

		String filename = "output.data";

		if (out == 0){
		Animate.main(args);
		} else if (out == 1){
		filename = "field.data";   
		PrintWriter output = new PrintWriter(new FileWriter(filename));
		Measurements.fieldMeasurements(c,p,W,output);
		} else if (out == 2){
		filename = "check.data";   
		PrintWriter output = new PrintWriter(new FileWriter(filename));
		Measurements.magCheck(c,p,W,output);
		} 

	}

}

	//----------------------------------------------------------------------

/**


		What's supposed to happen at the boundary all zero at the edge or beyond?
		***Think of a better boundary sln (might be in your notes)***
		Normalise the lengths of the vectors
		Check that the angle is in correct units
		am I okay initilising arrayPrime the way I do?


								  (1/6)*(array[i+1][j][k]+array[i==0?W:i-1][j][k]
									+array[i][j+1][k]+array[i][j==0?W:j-1][k] 
									+array[i][j][k+1]+array[i][j][k==0?W:k-1]
									+p[i][j][k]);


		Old jacobi		
					if(i==W-1 && j==W-1 ){arrayPrime[i][j] = 0.25*(array[i-1][j]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(i==W-1 && j==0){arrayPrime[i][j] = 0.25*(array[i-1][j]+array[i][j+1]+deltaX*deltaX*p[i][j]);
					}else if(i==0 && j==W-1){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(i==0 && j==0){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i][j+1]+deltaX*deltaX*p[i][j]);
					}else if(i==W-1){arrayPrime[i][j] = 0.25*(array[i-1][j]+array[i][j+1]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(j==W-1){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i-1][j]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else if(j==0){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i-1][j]+array[i][j+1]+deltaX*deltaX*p[i][j]);
					}else if(i==0){arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i][j+1]+array[i][j-1]+deltaX*deltaX*p[i][j]);
					}else{arrayPrime[i][j] = 0.25*(array[i+1][j]+array[i-1][j]+array[i][j+1]+array[i][j-1]+deltaX*deltaX*p[i][j]);}}}










**/	
