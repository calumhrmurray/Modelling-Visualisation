import java.awt.Color;
import java.awt.Frame;
import java.awt.Graphics2D;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;

class Animate {	// -------------------------------------------------------------

	static PDE Pmodel = new PDE();
	// k slice observing
	static int z = 25; // change this variable
	// number of pixels per arrow
	static int sq = 21;
	static int sq2 = 10;

	static final double tPI = 2.0*Math.PI;
	static void init(final BufferedImage bi, final int W) {
		Pmodel.update(W);
		double[][][] c = Pmodel.getPhiArray();
		int colour = 0;
		int xP = 0;
		int yP = 0;
		double maxPo = 0; 
		maxPo = maxPot(W);

		// arrows
		double maxLe = 0;
		maxLe = maxLen(W);
		final int w = bi.getWidth();
		final int h = bi.getHeight();
		final Graphics2D g = bi.createGraphics();
      		double[] field = new double[2];

		// print image	
		for (int x = 0; x < W; x++){
		for (int y = 0; y < W; y++){ colour = (int)Math.round(255*c[x][y][z]/maxPo); 
					     field = PDE.eField(c,x,y,z,W);
					     xP = x*sq;
					     yP = y*sq;		
					     //each point in the original array corresponds to a 21x21 drawn in the graphics
					     for (int i = 0; i < sq; i++){
					     for (int j = 0; j < sq; j++){ bi.setRGB(xP+i, yP+j, new Color(colour,colour,colour).getRGB());}} 
					     drawArrow(g, xP+sq2, yP+sq2, 5, field[1]); }}	
	}
  //5*field[0]/maxLe

	static void update(final BufferedImage bi, final int W) {
		Pmodel.jacobi(W); 
		double[][][] c = Pmodel.getPhiArray();
		int colour = 0;
		int xP = 0;
		int yP = 0;
		double maxPo = 0; 
		maxPo = maxPot(W);

		// arrows
		double maxLe = 0;
		maxLe = maxLen(W);
		final int w = bi.getWidth();
		final int h = bi.getHeight();
		final Graphics2D g = bi.createGraphics();
      		double[] field = new double[2];

		// print image	
		for (int x = 0; x < W; x++){
		for (int y = 0; y < W; y++){ colour = (int)Math.round(255*c[x][y][z]/maxPo); 
					     field = PDE.eField(c,x,y,z,W);
					     xP = x*sq;
					     yP = y*sq;		
					     for (int i = 0; i < sq; i++){
					     for (int j = 0; j < sq; j++){ bi.setRGB(xP+i, yP+j, new Color(colour,colour,colour).getRGB());}} 
					     drawArrow(g, xP, yP, 5, tPI); }}

	}

	//----------------------------------------------------------------------
	// scaling methods 

	// not currently neccesary to do over the entire array (just the slice we are concerned with)
	// leaving it general for now

	// find max val of potential
	static double maxPot(final int W) {
		double[][][] c = Pmodel.getPhiArray(); 
                double max = 0;               
                for(int i=0; i< W; i++){
		for(int j=0; j< W; j++){
		for(int k=0; k< W; k++){if(c[i][j][k] > max){max = c[i][j][k];}}}}
		return max;
	}

	// find max field vector values ignoring value at the point charge
	static double maxLen(final int W) {
		double[][][] c = Pmodel.getPhiArray(); 
		double maxL = 0; 
		double val = 0;
                for(int i=0; i< W; i++){
		for(int j=0; j< W; j++){
		for(int k=0; k< W; k++){ val = PDE.eField(c,i,j,k,W)[0];
					//if(i == (int)Math.round(W/2) & j == (int)Math.round(W/2) & k == (int)Math.round(W/2)){
					if(val > maxL){maxL = val;}
		}}}
		return maxL;
	}		

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 3) throw new Exception("Arguments: width[pixels] height[pixels] period[milliseconds]");
		final int W = Integer.parseInt(args[0]);
		final int acc = Integer.parseInt(args[1]);
		final int out = Integer.parseInt(args[2]);	
		
		final int imW = sq*W;
	
		final BufferedImage bi = new BufferedImage(imW, imW, BufferedImage.TYPE_INT_RGB);
		final Object lock = new Object();
		final Frame f = new Frame();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate " + args[0] + " " + args[1] + " " + args[2]);
		f.setVisible(true);
		f.setSize(imW, imW + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});

		init(bi,W);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {f.getGraphics().drawImage(bi, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}}, 0, 1);
		//new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {update(bi,W);}}}, 0, 1);
	}

 	   static final double qPI = 0.25*Math.PI;
 	   static void drawArrow(final Graphics2D g, final int x1, final int y1, final double l1, final double a) {
		g.setColor(Color.BLUE);
		final int x2 = x1 + (int)Math.round(l1*Math.cos(a));
		final int y2 = y1 + (int)Math.round(l1*Math.sin(a));
		g.drawLine(x1, y1, x2, y2);// tail
		final double l2 = 0.25*l1;
		final double am = a - qPI;
		final double ap = a + qPI;
		// tip:
		g.drawLine(x2, y2, x2 - (int)Math.round(l2*Math.cos(am)), y2 - (int)Math.round(l2*Math.sin(am)));
		g.drawLine(x2, y2, x2 - (int)Math.round(l2*Math.cos(ap)), y2 - (int)Math.round(l2*Math.sin(ap)));
 	   }
}











