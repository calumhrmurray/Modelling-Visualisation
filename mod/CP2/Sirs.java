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
		l = Rules.mcSweep(l,W,p1,p2,p3);
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
				        if (Math.random()<p){l[i][j]=1;}
				        if (Math.random()<p){l[i][j]=0;}
				        if (Math.random()<p){l[i][j]=-1;}}}
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
		if (args.length != 6) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final double p1 = Double.parseDouble(args[1]);
		final double p2 = Double.parseDouble(args[2]);
		final double p3 = Double.parseDouble(args[3]);
		final int out = Integer.parseInt(args[4]);
		final double im = Double.parseDouble(args[5]);

		l = new int[W][W];
		// Initialise model	
		initImmune(W,im);

		String filename = "output.data";

		// Choose output
		if ( out == 0 ){ 
		Animate.main(args);  // normal animation
		} else if ( out == 1){
		filename = "phaseDiagram.data";   
		PrintWriter output = new PrintWriter(new FileWriter(filename));
		Measurements.phaseDiagram(l,W,p2,output);	// phase diagram (model is iniated inside this model)
		} else if ( out == 2){
		filename = "infectedTimeDiagram.data";
		PrintWriter output = new PrintWriter(new FileWriter(filename));
		Measurements.infectedTimeDiagram(l,W,p1,p2,p3,output); // model is initiated here
		} else if ( out == 3){
		filename = "immunityDiagram.data";   
		PrintWriter output = new PrintWriter(new FileWriter(filename));
		Measurements.immunityDiagram(l,W,p2,p3,output);	// phase diagram (model is iniated inside this model)			
		}




	}

}


/** Check list

		check nearest neighbour calc (8 or 4?)
		check mcSweep
		\/\/\/\//\/not sure about the phase diagram - contour plot etc.\/\/\/\/\/
		\/\/\//\/\/\broken not changing susceptible ones to infected - issue with infection or nearestn\/\/\//\/
		\/\/\/\/\/check initialisation\/\/\//\/\/\/\
		why are p1 and p3 going up weirdly?
		\//\/\/\/\/include method for changeing number of immune cells into input\/\//\/\/\/\/\/
		-corresponding graph
		neaten up methods
		make sure all the graphs have error bars

*/

/** Ideas

		\/\/\/\/could draw the time vs infection graph in a graphics box as the simulation is running, would look cool\/\//\/\
		add start button
		switch between initial states before clicking start
		change % of immune cells before starting
		choose whether to draw graph or not
		add labels to graph legend	
	

*/











