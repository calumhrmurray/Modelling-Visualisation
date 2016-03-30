/*
 *	PDE solver
 */
import java.io.*;

public class PDE{

	private static double[][][] c;
	private static double[][][] p;
	private static double deltaX = 1;
	private static double perm = 1;
	private static double mu = 1;

	double[][][] getPhiArray() {
		return c;
	}

	double[][][] getDensityArray() {
		return p;
	}

	//----------------------------------------------------------------------
	// Charge distributions


	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] pointCharge(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}} 
		p[mid][mid][mid] = 1000000;
		return p;
	}

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] offCentrePointCharge(final int W){
		int mid = (int)Math.round(W/4);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}} 
		p[mid][mid][mid*2] = 10000;
		return p;
	}

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] quadrupole(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}}
		p[25][25][mid] = 10000;
		p[35][25][mid] = -10000;
		p[25][35][mid] = -10000;
		p[35][35][mid] = 10000;
		return p;
	}

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] gaussian(final int W){
		double sigma = .000001;
		double tpi = 2*Math.PI;
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = (1/(sigma*sigma*sigma*Math.sqrt(tpi*tpi*tpi)))*Math.exp((-1/(2*sigma*sigma))*((i-mid)*(i-mid)*(j-mid)*(j-mid)*(k-mid)*(k-mid)));}}}
		return p;
	}
	
	//----------------------------------------------------------------------
	// Current distributions

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] wire(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid && j==mid){
						p[i][j][k] = 10000; 
					} else {		
						p[i][j][k] = 0;
					}
					 }}} 
		return p;
	}

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] pwires(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid+10 && j==mid){
						p[i][j][k] = 10000; 
					} else if (i==mid-10 && j==mid){		
						p[i][j][k] = 10000; 
					} else {
						p[i][j][k] = 0;
					}
					 }}} 
		return p;
	}

	// adds an extra row column and row of zeros beyond the boundary
	static double[][][] apwires(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid+10 && j==mid){
						p[i][j][k] = 10000; 
					} else if (i==mid-10 && j==mid){		
						p[i][j][k] = -10000; 
					} else {
						p[i][j][k] = 0;
					}
					 }}} 
		return p;
	}
	

	//----------------------------------------------------------------------

	static double[][][] init(final int W){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){c[i][j][k] = 0;}}}
		return c;
	}

	static void update(final int W, final double acc, int choice){
		double error = acc+1;
		int i = 0;
		while(error>acc){
			i++;
			if (choice==0){
				error = jacobi(W);
				System.out.printf("Using: Jacobi  Trials: "+i+ "\r");
			} else if (choice==1){
				error = gaussSe(W);
				System.out.printf("Using: Gauss Seidel  Trials: "+i+" \r");
			} else if (choice==2){
				error = SORgauss(W,1.2);
				System.out.printf("Using: SOR  Trials: "+i+" \r");
			}
			if (i>100){
				break;
			}
		}	
	}

	//----------------------------------------------------------------------
	// algorithms 	
	
	// why didn't it work without the placemarker?
	static double jacobi(final int W){
		double[][][] arrayPrime = c;
		double placemarker = 0;
		double sum = 0;
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
		for (int k=0; k<W; k++){placemarker = (c[i+1][j][k]+c[i==0?W:i-1][j][k]   // if an array element on the boundary is called value will be zero
									+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
					sum += Math.abs(placemarker - arrayPrime[i][j][k]); 
					arrayPrime[i][j][k] = placemarker;
					;}}}
		c = arrayPrime;
		System.out.printf("Error: "+sum/(W*W*W)+" ");
		return sum/(W*W*W);
	}

	// is this correct for k?	
	static double gaussSe(final int W){
		double placemarker = 0;
		double sum = 0;

		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 0) placemarker=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								sum += Math.abs(placemarker-c[i][j][k]);
								c[i][j][k] = placemarker;
								}}} 
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 1) placemarker=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								sum += Math.abs(placemarker-c[i][j][k]);
								c[i][j][k] = placemarker;
								}}}
		
		System.out.printf("Error: "+sum/(W*W*W)+" ");
		return sum/(W*W*W);
	}


	static double SORgauss(final int W, double omega){
		double placemarker = 0;
		double sum = 0;

		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 0) placemarker=(1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k]
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								sum += Math.abs(placemarker-c[i][j][k]);
								c[i][j][k] = placemarker;
								}}} 
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 1) placemarker=(1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								sum += Math.abs(placemarker-c[i][j][k]);
								c[i][j][k] = placemarker;
								}}}

		System.out.printf("Error: "+sum/(W*W*W)+" ");
		return 2*sum/(W*W*W);		
	}

	//----------------------------------------------------------------------
	// Field calculators

	// takes in the array filled with potential outputs magnitude of the eField and the direction
	static double[] eField(double[][][] array, int i, int j, int k, int W){
		int mid = (int)Math.round(W/2);
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
		// think a bit more about the z*z terms
		if(y<0){
		field[1] = -Math.acos((x+z*z)/(Math.sqrt(x*x+y*y+z*z)*Math.sqrt(1+z*z)));
		} else if (y>0){
		field[1] = Math.acos((x+z*z)/(Math.sqrt(x*x+y*y+z*z)*Math.sqrt(1+z*z)));
		}
		return field;
		
	}

	// takes in the array filled with potential outputs magnitude of the mField and the direction
	static double[] mField(double[][][] array, int i, int j, int k, int W){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		double x = 0; double y = 0; double z = 0;
		// with boundary conditions
		x = (array[i][j+1][k]-array[i][j==0?W:j-1][k])/(2*deltaX);
		y = (array[i+1][j][k]-array[i==0?W:i-1][j][k])/(2*deltaX);
		// magnitude of the field
		// all need to be between 0 and 1
		// length	
		field[0] = Math.sqrt(x*x+y*y+z*z);
		// angle
		// think a bit more about the z*z terms
		if(y<0){
		field[1] = Math.acos((x+z*z)/(Math.sqrt(x*x+y*y+z*z)*Math.sqrt(1+z*z)));
		} else if (y>0){
		field[1] = -Math.acos((x+z*z)/(Math.sqrt(x*x+y*y+z*z)*Math.sqrt(1+z*z)));
		}
		return field;
		
	}

	//----------------------------------------------------------------------
	// Measurements

	static void fieldMeasurements(final int W, final double acc, int choice, PrintWriter output){
		int mid = (int)Math.round(W/2);
		update(W,acc,choice);	
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
	        output.printf(" "+mag(i,j,k,mid,mid,mid)+" "+c[i][j][k]+"\n");			
		}}
		output.printf("\n");}
	} 


	// |r| calculator
	static double mag(int i,int j,int k,int x,int y,int z){
		double r = 0;
		return Math.sqrt((i-x)*(i-x)+(j-y)*(j-y)+(k-z)*(k-z));
	}

	// find max field vector values, ignoring point charge
	static double maxLen(final int W) {
		int mid = (int)Math.round(W/2);
		double maxL = 0; 
		double val = 0;
                for(int i=0; i< W; i++){
		for(int j=0; j< W; j++){
		for(int k=0; k< W; k++){ val = PDE.eField(c,i,j,k,W)[0];
					if(i == mid & j == mid & k == mid){

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

		String filename = "output.data";

		if (out == 0){
		// Initialise model
		//pointCharge(W);
		//quadrupole(W);
		//gaussian(W);
		//offCentrePointCharge(W);
		//wire(W);
		//pwires(W);
		apwires(W);
		init(W);
		Animate.main(args);
		} else if (out == 1){
		// Initialise model
		pointCharge(W);
		init(W);
		filename = "field.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		fieldMeasurements(W,acc,0,output);
		} else if (out == 2){
		// Initialise model
		pointCharge(W);
		init(W);
		filename = "gfield.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		fieldMeasurements(W,acc,1,output);
		} else if (out == 3){
		// Initialise model
		quadrupole(W);
		init(W);
		filename = "qfield.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		fieldMeasurements(W,acc,1,output);
		} 
	}

}

	//----------------------------------------------------------------------



/**


		****What's supposed to happen at the boundary all zero at the edge or beyond?****
		***Think of a better boundary sln (might be in your notes)***
		***Normalise the lengths of the vectors****
		***Check that the angle is in correct units***
		****am I okay initilising arrayPrime the way I do?****
		sort out the for loop in the efield calculator
		****add acc thing*****
		*****sort out the angles!!****
		****move on to magnetics******
		*****laugh*******
		calculate convergence times
		get rid of arrrows at the point charges


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





	static void efieldMeasurements(final int W, PrintWriter output){
		for (int i=0;i<100;i++){
			PDE.jacobi(W);
		}
		double[] field = new double[2];
		int mid = (int)Math.round(W/2);
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		field = eField(c,i,j,mid,W);
	        output.printf(" "+mag(i,j,mid,mid,mid,mid)+" "+i+" "+j+" "+mid+" "+field[0]+" "+field[1]+"\n");			
		}}
	} 




**/	
