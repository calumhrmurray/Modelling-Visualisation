/**
 *  Measurements taken by PDE
 *
 */
import java.io.*;

public class Measurements{

	public static void fieldMeasurements(double[][][] c,double[][][] p,final int W, PrintWriter output){
		for (int i=0;i<1000;i++){
			PDE.jacobi(W);
		}		
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
	        output.printf(" "+mag(i,j,k,500,500,500)+" "+c[i][j][k]+"\n");			
		}}}
	} 

	// |r| calculator
	public static double mag(int i,int j,int k,int x,int y,int z){
		double r = 0;
		return Math.sqrt((i-x)*(i-x)+(j-y)*(j-y)+(k-z)*(k-z));
	}


	public static void magCheck(double[][][] c,double[][][] p,final int W, PrintWriter output){	
		// records |E| at |r| for every lattice site
		for (int i=0; i<W; i++){
		for (int j=0; j<W; j++){ 
		for (int k=0; k<W; k++){
	        output.printf(""+i+" "+j+" "+k+" "+mag(i,j,k,50,50,50)+"\n");			
		}}}
	} 

}
