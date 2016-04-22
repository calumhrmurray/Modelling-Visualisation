/*
 *	PDE solver
 */
import java.io.*;

public class PDE{

	private static double[][] f; // array stores fisher stuff	
	private static double[][] n;
	private static double deltaX = 1.; // distance corresponding to each discrete box
	private static double alpha = 1.;
	private static double D = 1.;
	private static double deltaT = 0.1;


	// getter method used to find the array c
	double[][] getFisherArray() {
		return f;
	}
	// getter method used to find the array c
	double[][] getnArray() {
		return n;
	}


	//----------------------------------------------------------------------
	// Droplet

	static double[][] droplet(final int W){
		int mid = (int)Math.round(W/2);
		// radius
		double R = 10;
		for (int i=0; i<W+1; i++){
		for (int j=0; j<W+1; j++){ 
					if(Math.sqrt((i-mid)*(i-mid)+(j-mid)*(j-mid))<R){f[i][j]=1;}else{f[i][j] = 0.0;}}} 
		return f;
	}

	static double[][] onedim(final int W){
		int Xo = (int)Math.round(W/10);
		for (int i=0; i<W+1; i++){
		for (int j=0; j<1; j++){
					if(i<Xo){n[i][j]=1.;}else{n[i][j] = 0.0;}}}
		return n;
	}

	static double[][] expOneDim(final int W, double k){
		for (int i=0; i<W+1; i++){
		for (int j=0; j<1; j++){n[i][j]=Math.exp(-k*i*deltaX);}}
		return n;
	}
 

	//----------------------------------------------------------------------

	static void FisherUpdate(final int W, final double acc, int alg){
		double error = acc+1;
		int i = 0;

		i++;
		if (alg==0){
			error = Fisher(W);
			System.out.printf("Using: Fisher Trials: "+i+ "  Error: "+error+" \r");
		} else if (alg==1){
			//error = gaussSe(W);
			System.out.printf("Using: Gauss Seidel  Trials: "+i+"  Error: "+error+" \r");
		}
		
		
		if (alg==0){
			System.out.printf("Using: Fisher  Trials: "+i+ "  Error: "+error+"\n");
		} else if (alg==1){
			System.out.printf("Using: Gauss Seidel  Trials: "+i+"  Error: "+error+" \n");
		} 
			
	}

	//----------------------------------------------------------------------
	// algorithms 	
	// each of these algorithms updates the array c and outputs an error value for termination of the update method



	static double Fisher(final int W){
		double val = 0; double valP  = 0; double sum = 0;

		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){
					  val = f[i][j];
					  valP= f[i][j]+deltaT*D*(f[i+1][j]+f[i==0?W:i-1][j]+f[i][j+1]+f[i][j==0?W:j-1]-4*f[i][j])/(deltaX*deltaX)
						       +deltaT*alpha*f[i][j]*(1.0-f[i][j]);
					 sum += Math.abs((val-valP));
					 f[i][j] = valP;
								}}
		
		double error = sum/(W*W);
		return error;
	}

	static double OneFisher(final int W){
		double val = 0; double valP  = 0; double sum = 0;
		double dx = 0;
		// starts at i=1 so that n[0][1] = 0 always
		for (int i=1; i<W; i++){
		for (int j=0; j<1; j++){if(i==W){dx=0;} else {dx=deltaT*D*(n[i+1][j]+n[i-1][j]-2*n[i][j])/(deltaX*deltaX);}
					 val = n[i][j];
					 valP= n[i][j]+dx+deltaT*alpha*n[i][j]*(1.0-n[i][j]);
					 sum += Math.abs((val-valP));
					 n[i][j] = valP;
								}}
		
		double error = sum/(W);
		return error;
	}
	
	//----------------------------------------------------------------------
	// Measurements

	static void Integral(final int W, PrintWriter output){
		double sum = 0;
		expOneDim(W,0.5);	
		// loops over a range of omega
		for(int t=0; t<1000; t++){
			sum = 0;
			// sum over array
			for (int i=1; i<W; i++){
			for (int j=0; j<1; j++){ sum += n[i][j];
					         }}
	       		output.printf(" "+t+" "+sum+"\n");
			OneFisher(W);
		}
		output.close();
	} 
	

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 5) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]); // size of system
		final double acc = Double.parseDouble(args[1]); // input for convergence
		final int out = Integer.parseInt(args[2]); // output choice	
		final int alg = Integer.parseInt(args[3]); // algorithm choice
		final int ch = Integer.parseInt(args[4]); // point charge/ quadrupole etc. choice
		
		f = new double[W+1][W+1];
		n = new double[W+1][1];

		String filename = "output.data";

		// initialises charge/current distribution depending on input ch
		if (ch==0){
		droplet(W); // initialises f
		} else if (ch==1){
		expOneDim(W,1);
		} 
		
		// initialises c
		//init(W);

		// decides output
		if (out == 0){
		// runs animation
		Animate.main(args);
		} else if (out == 1){
		// performs efield measurements 100 0.0001 1 2 0
		filename = "5Integral.data"; 
		PrintWriter output = new PrintWriter(new FileWriter(filename));  
		Integral(W,output);
		} 
	}

}

