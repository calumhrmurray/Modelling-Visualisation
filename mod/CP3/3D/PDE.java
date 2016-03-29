/*
 *	PDE solver
 */
import java.io.*;

public class PDE{

	private static double[][][] c;
	private static double[][][] p;
	private static double deltaX = 1;
	private static double perm = .1;

	double[][][] getPhiArray() {
		return c;
	}

	double[][][] getDensityArray() {
		return p;
	}

	//----------------------------------------------------------------------

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] pointCharge(final int W){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}}
		p[(int)Math.round(W/2)][(int)Math.round(W/2)][(int)Math.round(W/2)] = 10000;
		// why doesn't this appear
		p[40][40][40] = 10000;
		return p;
	}

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] quadrupole(final int W){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}}
		p[(int)Math.round(W/4)][(int)Math.round(W/4)][(int)Math.round(W/4)] = 10000;
		p[(int)Math.round(3*W/4)][(int)Math.round(W/4)][(int)Math.round(W/4)] = -10000;
		p[(int)Math.round(W/4)][(int)Math.round(3*W/4)][(int)Math.round(W/4)] = 10000;
		p[(int)Math.round(3*W/4)][(int)Math.round(3*W/4)][(int)Math.round(W/4)] = -10000;
		return p;
	}

	static double[][][] init(final int W){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){c[i][j][k] = 0;}}}
		return c;
	}

	static void update(final int W, final double acc){
		double error = acc;
					//System.out.printf("Error: "+error+"\n");	
		for (int i=0; i<2; i++){
					jacobi(W);
					//gaussSe(W);
					//System.out.printf("Error: "+error+"\n");	
					  }							
	}

	//----------------------------------------------------------------------

	static void jacobi(final int W){
		double[][][] arrayPrime = c;
		double sum = 0;
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
		for (int k=0; k<W; k++){arrayPrime[i][j][k] = (c[i+1][j][k]+c[i==0?W:i-1][j][k]   // if an array element on the boundary is called value will be zero
									+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
					;}}}
		sum = Math.sqrt((arrayPrime[26][26][26]-c[26][26][26])*(arrayPrime[26][26][26]-c[26][26][26]));
		c = arrayPrime;
		//System.out.printf("SUM:"+sum+"\n");
	}
	// is this correct for k?	
	static void gaussSe(final int W){
		// can this be done without arrayprime
		double[][][] arrayPrime = c;
		double sum = 0;
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 0) c[i][j][k]=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;}}} 
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 1) c[i][j][k]=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;}}}
		sum = Math.sqrt((arrayPrime[26][26][26]-c[26][26][26])*(arrayPrime[26][26][26]-c[26][26][26]));
		//System.out.printf("SUM:"+sum+"\n");
	}

	// takes in the array filled with potential outputs magnitude of the eField and the direction
	static double[] eField(double[][][] array, int i, int j, int k, int W){
      		double[] field = new double[2];
		double x = 0; double y = 0; double z = 0;
		// with boundary conditions
		x = -(array[i+1][j][k]-array[i==0?W:i-1][j][k])/(2*deltaX);
		y = -(array[i][j+1][k]-array[i][j==0?W:j-1][k])/(2*deltaX);
		z = -(array[i][j][k+1]-array[i][j][k==0?W:k-1])/(2*deltaX);
		// magnitude of the field
		// all need to be between 0 and 1
		// length	
		field[0] = Math.sqrt(x*x+y*y+z*z);
		// angle
		if(j<(int)Math.round(W/2)){
		field[1] = -Math.acos((x)/(Math.sqrt(x*x+y*y+z*z)));
		} else if (j>=(int)Math.round(W/2)){
		field[1] = Math.acos((x)/(Math.sqrt(x*x+y*y+z*z)));
		}
		return field;		
		
	}

	//----------------------------------------------------------------------
	// Measurements

	static void fieldMeasurements(final int W, PrintWriter output){
		for (int i=0;i<1000;i++){
			PDE.jacobi(W);
		}		
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
	        output.printf(" "+mag(i,j,k,(int)Math.round(W/2),(int)Math.round(W/2),(int)Math.round(W/2))+" "+c[i][j][k]+" "+i+" "+j+" "+k+"\n");			
		}}
		output.printf("\n");}
	} 

	static void hfieldMeasurements(final int W, PrintWriter output){
		for (int i=0;i<1000;i++){
			PDE.jacobi(W);
		}	
		
		int k = 50;
	
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
	        output.printf(" "+mag(i,j,k,(int)Math.round(W/2),(int)Math.round(W/2),(int)Math.round(W/2))+" "+c[i][j][k]+" "+i+" "+j+" "+k+"\n");			
		}
		output.printf("\n");}
	} 

	static void GafieldMeasurements(final int W, PrintWriter output){
		for (int i=0;i<1000;i++){
			PDE.gaussSe(W);
		}		
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
	        output.printf(" "+mag(i,j,k,(int)Math.round(W/2),(int)Math.round(W/2),(int)Math.round(W/2))+" "+c[i][j][k]+"\n");			
		}}
		output.printf("\n");}
	} 

	static void efieldMeasurements(final int W, PrintWriter output){
		for (int i=0;i<100;i++){
			PDE.jacobi(W);
		}
		double[] field = new double[2];
		int k = (int)Math.round(W/2);		
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		field = eField(c,i,j,k,W);
	        output.printf(" "+mag(i,j,k,(int)Math.round(W/2),(int)Math.round(W/2),(int)Math.round(W/2))+" "+i+" "+j+" "+k+" "+field[0]+" "+field[1]+"\n");			
		}}
	} 

	// |r| calculator
	static double mag(int i,int j,int k,int x,int y,int z){
		double r = 0;
		return Math.sqrt((i-x)*(i-x)+(j-y)*(j-y)+(k-z)*(k-z));
	}

	// find max field vector values
	static double maxLen(final int W) {
		double maxL = 0; 
		double val = 0;
                for(int i=0; i< W; i++){
		for(int j=0; j< W; j++){
		for(int k=0; k< W; k++){ val = PDE.eField(c,i,j,k,W)[0];
					if(i == (int)Math.round(W/2) & j == (int)Math.round(W/2) & k == (int)Math.round(W/2)){

					} else if(val > maxL){maxL = val;}
		}}}
		return maxL;
	}	

	static void maxLength(final int W	){
		for (int i=0;i<100;i++){
			PDE.jacobi(W);
		}
		double length = maxLen(W);
	        System.out.println("Max length:  "+length);		

	} 

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 3) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final double acc = Double.parseDouble(args[1]); // input for convergence
		final int out = Integer.parseInt(args[2]); // output choice

		c = new double[W+1][W+1][W+1];
		p = new double[W+1][W+1][W+1];		

		// Initialise model
		//pointCharge(W);
		quadrupole(W);
		init(W);

		String filename = "output.data";

		if (out == 0){
		Animate.main(args);
		} else if (out == 1){
		filename = "field.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		fieldMeasurements(W,output);
		} else if (out == 2){
		filename = "vector.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		efieldMeasurements(W,output);
		} else if (out == 3){
		maxLength(W);
		} else if (out == 4){
		filename = "Gfield.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		GafieldMeasurements(W,output);
		} else if (out == 5){
		filename = "hfield.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		hfieldMeasurements(W,output);
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
		sort out the for loop in the efield calculator


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
