/*
 *	PDE solver
 */
import java.io.*;

public class PDE{

	private static double[][][] c; // array stores the potential
	private static double[][][] p; // arrray stores charge or current density
	private static double deltaX = 1; // distance corresponding to each discrete box
	private static double perm = 1; // dielectric constant of the vacuum
	private static double mu = 1; // magnetic constant of permeability

	// getter method used to find the array c
	double[][][] getPhiArray() {
		return c;
	}
	// getter method used to find the array p
	double[][][] getDensityArray() {
		return p;
	}

	//----------------------------------------------------------------------
	// Charge distributions
	// These methods initiate each of the charge distributions - just numbers in an array	
	// an extra row column and row of zeros beyond the boundary to be called in the algorithms	

	static double[][][] pointCharge(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0.0;}}} 
		p[mid][mid][mid] = 10000.0;
		return p;
	}


	static double[][][] offCentrePointCharge(final int W){
		int quarter = (int)Math.round(W/4);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0.0;}}} 
		p[quarter][quarter][quarter*2] = 10000.0;
		return p;
	}


	static double[][][] quadrupole(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = 0;}}}
		p[mid+5][mid][mid] = -10000.0;
		p[mid-5][mid][mid] = 10000.0;
		p[mid-5][mid+5][mid] = -10000.0;
		p[mid+5][mid+5][mid] = 10000.0;
		return p;
	}


	static double[][][] gaussian(final int W){
		double sigma = .0000000000001;
		double tpi = 2*Math.PI;
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){p[i][j][k] = (1/(sigma*sigma*sigma*Math.sqrt(tpi*tpi*tpi)))*Math.exp((-1/(2*sigma*sigma))*((i-mid)*(i-mid)*(j-mid)*(j-mid)*(k-mid)*(k-mid)));}}}
		return p;
	}
	
	//----------------------------------------------------------------------
	// Current distributions
	// These methods initiate each of the current distributions - just numbers in an array	
	// an extra row column and row of zeros beyond the boundary is added to be called in the algorithms	

	static double[][][] wire(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid && j==mid){
						p[i][j][k] = 10000.0; 
					} else {		
						p[i][j][k] = 0.0;
					}
					 }}} 
		return p;
	}


	static double[][][] pwires(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid+10 && j==mid){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-10 && j==mid){		
						p[i][j][k] = 10000.0; 
					} else {
						p[i][j][k] = 0.0;
					}
					 }}} 
		return p;
	}


	static double[][][] apwires(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid+10 && j==mid){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-10 && j==mid){		
						p[i][j][k] = -10000.0; 
					} else {
						p[i][j][k] = 0.0;
					}
					 }}} 
		return p;
	}

	static double[][][] mpwires(final int W){
		int mid = (int)Math.round(W/2);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){ if (i==mid+5 && j==mid){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-5 && j==mid){		
						p[i][j][k] = -10000.0; 
					} else if (i==mid+5 && j==mid+10){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-5 && j==mid+10){		
						p[i][j][k] = -10000.0; 
					} else if (i==mid+5 && j==mid-10){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-5 && j==mid-10){		
						p[i][j][k] = -10000.0; 
					} else if (i==mid+5 && j==mid+20){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-5 && j==mid+20){		
						p[i][j][k] = -10000.0; 
					} else if (i==mid+5 && j==mid-20){
						p[i][j][k] = 10000.0; 
					} else if (i==mid-5 && j==mid-20){		
						p[i][j][k] = -10000.0; 
					} else {
						p[i][j][k] = 0.0;
					}
					 }}} 
		return p;
	}
	

	//----------------------------------------------------------------------

	// initiates the electric potential/ magnetic potential array
	static double[][][] init(final int W){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
		for (int k=0; k<W+1; k++){c[i][j][k] = 0;}}}
		return c;
	}

	// update method
	// used in all instances to run the algorithms
	// alg = choice of algorithm
	// terminates when error<acc or system is deemed divergent (i> large value) 
	// this method also prints to the terminal to confirm which algorithms is being used
	static void update(final int W, final double acc, int alg){
		double error = acc+1;
		int i = 0;
		while(error>acc){
			i++;
			if (alg==0){
				error = jacobi(W);
				System.out.printf("Using: Jacobi  Trials: "+i+ "  Error: "+error+" \r");
			} else if (alg==1){
				error = gaussSe(W);
				System.out.printf("Using: Gauss Seidel  Trials: "+i+"  Error: "+error+" \r");
			} else if (alg==2){
				error = SORgauss(W,1.4);
				System.out.printf("Using: SOR  Trials: "+i+"  Error: "+error+" \r");
			}
			if (i>10000){
				break;
			}
		}
		
		if (alg==0){
			System.out.printf("Using: Jacobi  Trials: "+i+ "  Error: "+error+"\n");
		} else if (alg==1){
			System.out.printf("Using: Gauss Seidel  Trials: "+i+"  Error: "+error+" \n");
		} else if (alg==2){
			System.out.printf("Using: SOR  Trials: "+i+"  Error: "+error+"\n");
		}
			
	}

	//----------------------------------------------------------------------
	// algorithms 	
	// each of these algorithms updates the array c and outputs an error value for termination of the update method
	
	//Jacobi algorithm
	static double jacobi(final int W){
		double[][][] arrayPrime = c;
		double val = 0; double valP  = 0; double sum = 0;
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
		for (int k=0; k<W; k++){val = c[i][j][k];
					valP = (c[i+1][j][k]+c[i==0?W:i-1][j][k]   // if an array element on the boundary is called value will be zero
									+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
					sum += Math.abs((val-valP));
					//sum += Math.abs(c[i][j][k]-arrayPrime[i][j][k]);
					arrayPrime[i][j][k] = valP;
					}}}
		c = arrayPrime;
		double error = sum/(W*W*W);
		return error;
	}

	//Gauss Seidel algorithm
	static double gaussSe(final int W){
		double val = 0; double valP  = 0; double sum = 0;

		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){ val = c[i][j][k];
					  valP=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
						+c[i][j][k+1]+c[i][j][k==0?W:k-1]
						+deltaX*deltaX*p[i][j][k])/6;
					 sum += Math.abs((val-valP));
					 c[i][j][k] = valP;
								}}}
		
		double error = sum/(W*W*W);
		return error;
	}

	// SOR algorithm
	static double SORgauss(final int W, double omega){
		double val = 0; double valP  = 0; double sum = 0;

		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){ val = c[i][j][k];
					  valP=(1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
						+c[i][j][k+1]+c[i][j][k==0?W:k-1]
						+deltaX*deltaX*p[i][j][k])/6;
					 sum += Math.abs((val-valP));
					 c[i][j][k] = valP;
								}}}


		double error = sum/(W*W*W);
		return error;
	}

	//----------------------------------------------------------------------
	// Field calculators

	// Electric Field calculator
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
		// this is required for the measurements but would be incorrect to use for this visualisation as it has the z term
		field[0] = Math.sqrt(x*x+y*y+z*z);
		// angle
		if(y<0){
		field[1] = -Math.acos((x)/(Math.sqrt(x*x+y*y)));
		} else if (y>0){
		field[1] = Math.acos((x)/(Math.sqrt(x*x+y*y) ));
		}
		return field;		
	}

	// Magnetic Field calculator
	// takes in the array filled with potential outputs magnitude of the mField and the direction
	static double[] mField(double[][][] array, int i, int j, int k, int W){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		double x = 0; double y = 0; double z = 0;
		// with boundary conditions
		x = (array[i][j+1][k]-array[i][j==0?W:j-1][k])/(2*deltaX);
		y = (array[i+1][j][k]-array[i==0?W:i-1][j][k])/(2*deltaX);
		// magnitude of the field
		field[0] = Math.sqrt(x*x+y*y);
		// angle
		if(y<0){
		//field[1] = Math.acos((x+z*z)/(Math.sqrt(x*x+y*y+z*z)*Math.sqrt(1+z*z)));
		field[1] = Math.acos((x)/(Math.sqrt(x*x+y*y)));
		} else if (y>0){
		//field[1] = -Math.acos((x+z*z)/(Math.sqrt(x*x+y*y+z*z)*Math.sqrt(1+z*z)));
		field[1] = -Math.acos((x)/(Math.sqrt(x*x+y*y+z*z)));
		}
		return field;
		
	}

	//----------------------------------------------------------------------
	// Measurements

	// Ouputs the electric field and potential at radius |r|
	static void efieldMeasurements(final int W, final double acc, final int alg, PrintWriter output){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		update(W,acc,alg);	
		// records potential at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
					     field = PDE.eField(c,i,j,k,W);
	       				     output.printf(" "+mag(i,j,k,mid,mid,mid)+" "+field[0]+" "+c[i][j][k]+"\n");			
		}}
		output.printf("\n");}
	} 

	// Ouputs the electric field and potential at radius |r|
	static void refieldMeasurements(final int W, final double acc, final int alg, PrintWriter output){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		update(W,acc,alg);	
		// records potential at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
					     field = PDE.eField(c,i,j,k,W);
						if (i<50){
	       				     output.printf(" "+-mag(i,j,k,mid,mid,mid)+" "+field[0]+" "+c[i][j][k]+"\n");
						} else {
	       				     output.printf(" "+mag(i,j,k,mid,mid,mid)+" "+field[0]+" "+c[i][j][k]+"\n");
						}			
		}}
		output.printf("\n");}
	} 

	// Ouputs the magnetic field at radius |r|
	static void mfieldMeasurements(final int W, final double acc, final int alg, PrintWriter output){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		update(W,acc,alg);	
		// records potential at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
					     field = PDE.mField(c,i,j,mid,W);	
	       				     output.printf(" "+mag(i,j,mid,mid,mid,mid)+" "+field[0]+" "+c[i][j][mid]+"\n");			
		}
		output.printf("\n");}
	} 

	// Ouputs the magnetic field at radius |r|
	static void solenoid(final int W, final double acc, final int alg, PrintWriter output){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		update(W,acc,alg);	
		// records potential at |r| for every lattice site
		for (int k=0; k<W; k++){
					     field = PDE.mField(c,mid,mid,k,W);	
	       				     output.printf(" "+mag(mid,mid,k,mid,mid,mid)+" "+field[0]+" "+c[mid][mid][k]+"\n");			
		}
		output.close();


	} 

	// Ouputs the electric field and potential at k=mid contour plot
	static void ContourEfield(final int W, final double acc, final int alg, PrintWriter output){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		update(W,acc,alg);	
		// records potential at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
					     field = PDE.eField(c,i,j,mid,W);
	       				     output.printf(" "+i+" "+j+" "+field[0]+" "+c[i][j][mid]+"\n");			
		}
		output.printf("\n");}
	} 

	// Ouputs the electric field and potential at k=mid contour plot
	static void ContourMfield(final int W, final double acc, final int alg, PrintWriter output){
		int mid = (int)Math.round(W/2);
      		double[] field = new double[2];
		update(W,acc,alg);	
		// records potential at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
					     field = PDE.mField(c,i,j,mid,W);
	       				     output.printf(" "+i+" "+j+" "+field[0]+" "+c[i][j][mid]+"\n");			
		}
		output.printf("\n");}
	} 

	// SOR convergence test
	static void SORconvergence(final int W, final double acc, PrintWriter output){
		double omega = 1.9;
		double error = acc+1;	
		int counter = 0;
		// loops over a range of omega
		for(int i=0; i<20; i++){
			pointCharge(W);		
			init(W);
			counter = 0;	
			// same as with the update loop
			while(error>acc){
				counter++;
				error = SORgauss(W,omega);
				System.out.printf("Using: SOR  Omega: "+omega+"  Counter: "+counter+"\r");	
				// consider values that do not converge after 10000 runs to be divergent
				if (i>10000){
					// want to return NaN to show none convergence
					break;
				}
			}
			System.out.printf("\n");
	       		output.printf(" "+omega+" "+counter+" "+error+"\n");
			error = acc+1;	
			omega += 0.005;	
		}
		output.close();
	} 


	// method used to run update whilst checking that errors seem to be working correctly (observing the output to the terminal)
	static void conRun(final int W, final double acc, int alg){
		int mid = (int)Math.round(W/2);
		update(W,acc,alg);	

	}

	//----------------------------------------------------------------------
	// Calculators

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

	//find max field vector values
	static void maxLength(final int W){
		for (int i=0;i<100;i++){
			PDE.jacobi(W);
		}
		double length = maxLen(W);
	        System.out.println("Max length:  "+length);		

	} 

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 5) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]); // size of system
		final double acc = Double.parseDouble(args[1]); // input for convergence
		final int out = Integer.parseInt(args[2]); // output choice	
		final int alg = Integer.parseInt(args[3]); // algorithm choice
		final int ch = Integer.parseInt(args[4]); // point charge/ quadrupole etc. choice
		
		c = new double[W+1][W+1][W+1];
		p = new double[W+1][W+1][W+1];		

		String filename = "output.data";

		// initialises charge/current distribution depending on input ch
		if (ch==0){
		pointCharge(W);
		} else if (ch==1){
		quadrupole(W);
		} else if (ch==2){
		gaussian(W);
		} else if (ch==3){
		offCentrePointCharge(W);
		} else if (ch==4){
		wire(W);
		} else if (ch==5){
		pwires(W);
		} else if (ch==6){
		apwires(W);
		} else if (ch==7){
		mpwires(W);
		}
		
		// initialises c
		init(W);

		// decides output
		if (out == 0){
		// runs animation
		Animate.main(args);
		} else if (out == 1){
		// performs efield measurements 100 0.0001 1 2 0
		filename = "PointChargeE.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		efieldMeasurements(W,acc,alg,output);
		} else if (out == 2){
		// performs efield measurements on gaussian 100 0.0001 2 2 2
		filename = "Gaussian.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		efieldMeasurements(W,acc,alg,output);
		}else if (out == 3){
		// performs mfield measurements on slice 100 0.0001 3 2 7
		filename = "SolenoidCentral.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		solenoid(W,acc,alg,output);
		} else if (out == 4){
		// convergence of SOR for varying omega
		filename = "SORconvergence.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		SORconvergence(W,acc,output);
		} else if (out == 5){
		conRun(W,acc,alg);
		} else if (out == 6){
		// performs efield measurements 100 0.0001 1 2 0
		filename = "quadrupoleR.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		efieldMeasurements(W,acc,alg,output);
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
		****sort out the for loop in the efield calculator*****
		****add acc thing*****
		*****sort out the angles!!****
		****move on to magnetics******
		*****laugh*******
		******calculate convergence times******
		get rid of arrrows at the point charges
		******difficulty assigning double array to a value******
		*****check in the quadrupole that field lines point in the correct direction*****


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



		/*
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 0) c[i][j][k]=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								}}} 
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 1) c[i][j][k]=(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								}}}*/

		/*
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
		for (int k=0; k<W; k++){if ((i+j+k) % 2 == 0)  val = c[i][j][k];
							       valP = (1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k]
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
						sum += Math.abs((val-valP));
						c[i][j][k] = valP;
								}}} 

		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
		for (int k=0; k<W; k++){if ((i+j+k) % 2 == 1)  val = c[i][j][k];
								valP = (1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k]
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
						sum += Math.abs((val-valP));
						c[i][j][k] = valP;
								}}} */

		/*for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 0) c[i][j][k]=(1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k]
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								}}} 
		for (int i=0; i<W; i+=1){
		for (int j=0; j<W; j+=1){
		for (int k=0; k<W; k+=1){if ((i+j+k) % 2 != 1) c[i][j][k]=(1-omega)*c[i][j][k]+omega*(c[i+1][j][k]+c[i==0?W:i-1][j][k]+c[i][j+1][k]+c[i][j==0?W:j-1][k] 
									+c[i][j][k+1]+c[i][j][k==0?W:k-1]
									+deltaX*deltaX*p[i][j][k])/6;
								}}}*/
