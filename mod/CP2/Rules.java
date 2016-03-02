/**
 *	Rules, Static methods used to determine game rules
 */
import java.lang.Math;

public class Rules{

	//----------------------------------------------------------------------

	public static int[][] parallelUpdate(int[][] l, int W, double p1, double p2, double p3){
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){l = Life(l,W,i,j,p1,p2,p3);}}
	return l;
	}

	public static int[][] mcSweep(int[][] l, int W, double p1, double p2, double p3){
		// randomly choose N lattice sites update lattice after each time
		for(int x=0; x<W*W; x++){
					int i = (int)(Math.random()*W);
					int j = (int)(Math.random()*W);
					// is this correct?
					l = Life(l,W,i,j,p1,p2,p3);
					}
	
	return l;	
	}

	// Rules of SIRS
	public static int[][] Life(int[][] l, int W, int i, int j, double p1, double p2, double p3){
			int a = l[i][j];
			if (Infection(l,i,j,W)!=0 && a==1){
				if (Math.random() < p1){ l[i][j] = -1;}			// infected
			} else if (a ==-1){
				if (Math.random() < p2){ l[i][j] =  0;}			// recovers
			} else if (a== 0){
				if (Math.random() < p3){ l[i][j] =  1;}			// becomes susceptible
			} 
		return l;
	}
	
	//used to correct BCs
	static int mod (int x, int m){
		//m = Math.abs(m);
		return ( x % m + m) % m;
	}

	// if int n > 0 then a site has an infected nearest neighbour
	// there is no need to check every nearest neighbour if one infected neighbour is found
	// correct nearest neighbours or back to Ising nearest neighbours?
	static int Infection(int[][] l, int x, int y, int W){
		int n = 0;
		// mod corrects for boundary conditions
		if (l[x==W-1?0:(x+1)][y]==-1){
		n++;		
		}
		else if (l[x][y==W-1?0:(y+1)]==-1){
		n++;		
		}
		else if (l[x==0?W-1:(x-1)][y]==-1){
		n++;		
		}
		else if (l[x][y==0?W-1:(y-1)]==-1){
		n++;		
		}

		return n;	
	}

	// if int n > 0 then a site has an infected nearest neighbour
	// there is no need to check every nearest neighbour if one infected neighbour is found
	// correct nearest neighbours or back to Ising nearest neighbours?
	static int NearestNeighbour(int[][] l, int x, int y, int W){
		int n = 0;
		// mod corrects for boundary conditions
		if (l[mod(x+1,W)][y]==-1){
		n++;		
		}
		else if (l[mod(x+1,W)][mod(y+1,W)]==-1){
		n++;		
		}
		else if (l[x][mod(y+1,W)]==-1){
		n++;		
		}
		else if (l[mod(x-1,W)][mod(y+1,W)]==-1){
		n++;		
		}
		else if (l[mod(x-1,W)][y]==-1){
		n++;		
		}
		else if (l[mod(x-1,W)][mod(y-1,W)]==-1){
		n++;		
		}
		else if (l[x][mod(y-1,W)]==-1){
		n++;		
		}
		else if (l[mod(x+1,W)][mod(y-1,W)]==-1){
		n++;		
		}

		return n;	
	}

	//----------------------------------------------------------------------

	// Rules of SIRS
	public static int[][] ParalellLife(int[][] l, int W, double p1, double p2, double p3){
		// Apply the rules of Life
		for (int i = 0; i < W; i++){
		for (int j = 0; j < W; j++){
			if (Infection(l,i,j,W)>0 && l[i][j]==1){
				if (Math.random() < p1){ l[i][j] = -1;} 			// infected
			} else if (l[i][j]==-1){
				if (Math.random() < p2){ l[i][j] =  0;}			        // recovers
			} else if (l[i][j]==0){
				if (Math.random() < p3){ l[i][j] =  1;}	            		// becomes susceptible
			} 
		}}

		return l;
	}
	


}
