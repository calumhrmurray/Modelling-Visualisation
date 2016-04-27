/** runs a bash script */

import java.io.*;
import java.util.*;

class dynamicplot{

	public static void startPlot(String x, String y1, String y2, String y3) throws Exception{
		Process p = new ProcessBuilder("./start_dynamicplot.sh",x,y1,y2,y3).start();
	}		

	public static void updatePlot(String x, String y1, String y2, String y3) throws Exception {	
		Process p = new ProcessBuilder("./dynamicplot.sh",x,y1,y2,y3).start();
    	}

	public static void resetPlot() throws Exception {	
		Process p = new ProcessBuilder("./reset_dynamicplot.sh").start();
    	}

	public static void exitPlot() throws Exception {	
		Process p = new ProcessBuilder("./exit_dynamicplot.sh").start();
    	}


}
