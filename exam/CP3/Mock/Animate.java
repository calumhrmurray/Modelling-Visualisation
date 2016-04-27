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

	// creates buffered image displaying potential and field lines
	static void init(final BufferedImage bi, final int W, final double acc, final int alg, final int ch) {
		// initialises model
		Pmodel.droplet(W);
		// runs update
		Pmodel.FisherUpdate(W,acc,alg);
		// gets c the potential array
		double[][] f = Pmodel.getFisherArray();
		int colour = 0;
		double maxPo = 0;
		// finds potential maximum used for the colour scale 
		maxPo = maxPot(W);
		double minPo = 0;
		// finds potential minimum used for the colour scale 
		minPo = minPot(W);

		// print image	
		for (int x = 0; x < W; x++){
		for (int y = 0; y < W; y++){ // colour of a square
					     colour = (int)Math.round(255*(f[x][y]-minPo)/(maxPo-minPo)); 
					     bi.setRGB(x, y, new Color(colour,colour,colour).getRGB());}} 
	}

	static void update(final BufferedImage bi, final int W, final double acc, final int alg, final int ch) {		
		Pmodel.FisherUpdate(W,acc,alg);
		// gets c the potential array
		double[][] f = Pmodel.getFisherArray();
		int colour = 0;
		double maxPo = 0;
		// finds potential maximum used for the colour scale 
		maxPo = maxPot(W);
		double minPo = 0;
		// finds potential minimum used for the colour scale 
		minPo = minPot(W);

		// print image	
		for (int x = 0; x < W; x++){
		for (int y = 0; y < W; y++){ // colour of a square
					     colour = (int)Math.round(255*(f[x][y]-minPo)/(maxPo-minPo)); 
					     bi.setRGB(x, y, new Color(colour,colour,colour).getRGB());}} 
	}



	
	//----------------------------------------------------------------------
	// scaling methods 

	// not currently neccesary to do over the entire array (just the slice we are concerned with)
	// leaving it general for now

	// find max val of potential
	static double maxPot(final int W) {
		double[][] f = Pmodel.getFisherArray(); 
                double max = 0;               
                for(int i=0; i< W; i++){
		for(int j=0; j< W; j++){
					if(f[i][j] > max){max = f[i][j];}}}
		return max;
	}

	static double minPot(final int W) {
		double[][] f = Pmodel.getFisherArray(); 
                double min = 0;               
                for(int i=0; i< W; i++){
		for(int j=0; j< W; j++){
					if(f[i][j] < min){min = f[i][j];}}}
		return min;
	}
	

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 5) throw new Exception("Arguments: width[pixels] height[pixels] period[milliseconds]");
		final int W = Integer.parseInt(args[0]);
		final double acc = Double.parseDouble(args[1]); // input for convergence
		final int out = Integer.parseInt(args[2]); // output choice	
		final int alg = Integer.parseInt(args[3]); // algorithm choice
		final int ch = Integer.parseInt(args[4]); // point charge/ quadrupole etc. choice
		

		final BufferedImage bi = new BufferedImage(W, W, BufferedImage.TYPE_INT_RGB);
		final Object lock = new Object();
		final Frame f = new Frame();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate " + args[0] + " " + args[1] + " " + args[2]);
		f.setVisible(true);
		f.setSize(W, W + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});

		init(bi,W,acc,alg,ch);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {f.getGraphics().drawImage(bi, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}}, 0, 1);
		// no need to run updates as model is only display once error is below input error level
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {update(bi,W,acc,alg,ch);}}}, 0, 100);
	}

}
