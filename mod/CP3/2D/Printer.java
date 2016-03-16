/*
 *	SIRS model
 */
import java.io.*;

public class Printer{

	static PDE Pmodel = new PDE();

	static void init(final int W) {
		double[][] c = Pmodel.init(W);
		for (int i = 0; i < W; i++){for (int j = 0; j < W; j++){ System.out.printf(" "+(int)(c[i][j]+0.5d)+"     ");} System.out.printf("\n");}		
		System.out.printf("\n");	
	}

	static void update(final int W) {
		double[][] c = Pmodel.update(W); 
		for (int i = 0; i < W; i++){for (int j = 0; j < W; j++){ System.out.printf(" "+(int)(c[i][j]+0.5d)+"     ");} System.out.printf("\n");}		
		System.out.printf("\n");
	}

	static void updateB(final int W) {
		double[][] c = Pmodel.update(W); 
	}

	public static void main(final String[] args) throws Exception {
		if (args.length != 3) throw new Exception("Arguments: width[pixels] height[pixels] period[milliseconds]");
		final int W = Integer.parseInt(args[0]);
		final int acc = Integer.parseInt(args[1]);
		final int out = Integer.parseInt(args[2]);	
	
		init(W);
		for(int i = 0; i < 1000000;i++){updateB(W);}

		update(W);
		System.out.printf("Fin! \n");
	}

}
