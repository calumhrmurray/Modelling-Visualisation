/*
 *	SIRS model
 */
import java.io.*;

public class Sirs{

	private static int[][] l;

	int[][] getArray() {
		return l;
	}

	//----------------------------------------------------------------------

	static int[][] update(final int W, final double p1, final double p2, final double p3){
		return Rules.mcSweep(l,W,p1,p2,p3);
		// different update methods can be placed here such as parralel 
	}

	static int[][] updateAnimate(final int W, final double p1, final double p2, final double p3){
		int i = (int)(Math.random()*W);
		int j = (int)(Math.random()*W);
		l = Rules.Life(l,W,i,j,p1,p2,p3);
	return l;
	}

	
	// randomly initialises lattice equally S,I,R	
	static int[][] init(final int W) {
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ l[i][j] = (int) (Math.random()*3) - 1; }}
			
         return l;	
	}

	// adds Immune states = 2 with probability Im
	static int[][] initImmune(final int W, double Im){
		double p = (1-Im)/3;
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ l[i][j] = 0; // default susceptible (do not want unitialised state)
					if (Math.random()<Im){l[i][j]=2;} // immune
					else {l[i][j] = (int) (Math.random()*3) - 1;}
					}}
	return l;
	}

	// adds one infected square the rest are susceptible
	static int[][] initSus(final int W){
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ if (i==50 && j==50){l[i][j]=-1;} else {l[i][j]=1;}}} 
	return l;
	}

	//----------------------------------------------------------------------


	public static void main(final String[] args) throws Exception {
		if (args.length != 4) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final double p1 = Double.parseDouble(args[1]);
		final double p2 = Double.parseDouble(args[2]);
		final double p3 = Double.parseDouble(args[3]);

		l = new int[W][W];
		// Initialise model	
		init(W);

		String filename = "output.data";

		//Animate.main(args);
		//Animate1.main(args);
		Animate2.main(args);
		
		

	}

}












