/*
 *
 */
import java.io.*;

public class Game{

	private static int[][] l;

	int[][] getArray() {
		return l;
	}

	//----------------------------------------------------------------------

	static int[][] update(int W){
		l = Rules.Life(l,W);
	return l;
	}
	static int[][] init(final int W, final int c) {
		if (c == 0){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){l[i][j] = 1;}}
        	} else if (c==1){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){l[i][j] = -1;}}
		} else if (c==2){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){if (i==49 && j==50 ||i==50 && j==50 ||i==51 && j==50 ){l[i][j] = 1;} else{ l[i][j] = -1;}}}
		} else if (c==3){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){if (i==49 && j==50 ||i==49 && j==51 ||i==50 && j==50 ||i==50 && j==51){l[i][j] = 1;} else{ l[i][j] = -1;}}}
		} else if (c==4){
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){if (i==49 && j==50 ||i==50 && j==50 ||i==51 && j==50 ||i==51 && j==51||i==50 && j==52 
							){l[i][j] = 1;} else{ l[i][j] = -1;}}}
		} else {
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){if (Math.random() < 0.2){ l[i][j] = 1;} else{ l[i][j] = -1;}}}}
         return l;	
	}

	//----------------------------------------------------------------------


	public static void main(final String[] args) throws Exception {
		if (args.length != 2) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]);
		final int c = Integer.parseInt	(args[1]);


		String filename = "random.data";
		// start printer
		PrintWriter output = new PrintWriter(new FileWriter(filename));

		l = new int[W][W];
		// Initialise model	
		init(W,c);

		Animate.main(args);
	

		}




}