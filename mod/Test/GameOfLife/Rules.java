/**
 *	Rules, Static methods used to determine game rules
 */
import java.lang.Math;

public class Rules{

	// Rules of the game of life, for some reason not right probably an issue with getNcount
	public static int[][] Life(int[][] l, int W){
		// Apply the rules of Life
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
			// dies
			if (getNcount(l,i,j,W)<2 || getNcount(l,i,j,W)>3){
				l[i][j] = -1;
			}
			// born or stays alive
			else if (getNcount(l,i,j,W)==3){
				l[i][j] = 1;
			} 
		}}

		return l;
	}
	
	//used to correct BCs
	public static int mod (int x, int m){
		m = Math.abs(m);
		return ( x % m + m) % m;
	}

	public static int getNcount(int[][] l, int x, int y, int W){
		int n = 0;
		// mod corrects for boundary conditions
		if (l[mod(x+1,W)][y]==1){
		n++;		
		}
		if (l[mod(x+1,W)][mod(y+1,W)]==1){
		n++;		
		}
		if (l[x][mod(y+1,W)]==1){
		n++;		
		}
		if (l[mod(x-1,W)][mod(y+1,W)]==1){
		n++;		
		}
		if (l[mod(x-1,W)][y]==1){
		n++;		
		}
		if (l[mod(x-1,W)][mod(y-1,W)]==1){
		n++;		
		}
		if (l[x][mod(y-1,W)]==1){
		n++;		
		}
		if (l[mod(x+1,W)][mod(y-1,W)]==1){
		n++;		
		}

		return n;	
	}
	


}
