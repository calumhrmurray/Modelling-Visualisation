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
	static int z = 25;
	// number of pixels per arrow
	static int sq = 11;

	static void init(final BufferedImage bi, final int W) {
		double[][][] c = Pmodel.init(W);
		int colour = 0;
		int xP = 0;
		int yP = 0;
		// print image	
		for (int x = 0; x < W; x++){
		for (int y = 0; y < W; y++){ colour = (int)(c[x][y][z]+0.5d); 
					     xP = x*sq;
					     yP = y*sq;		
					     for (int i = 0; i < sq; i++){
					     for (int j = 0; j < sq; j++){ bi.setRGB(xP+i, yP+j, new Color(0,0,colour%255).getRGB());}} }}
	}

	static final double tPI = 2.0*Math.PI;
	static void update(final BufferedImage bi, final int W) {
		double[][][] c = Pmodel.update(W); 
		int colour = 0;
		int xP = 0;
		int yP = 0;
		// print image	
		for (int x = 0; x < W; x++){
		for (int y = 0; y < W; y++){ colour = (int)(c[x][y][z]+0.5d); 
					     xP = x*sq;
					     yP = y*sq;		
					     for (int i = 0; i < sq; i++){
					     for (int j = 0; j < sq; j++){ bi.setRGB(xP+i, yP+j, new Color(0,0,colour%255).getRGB());}} }}
		
		// arrows
		/*final int w = bi.getWidth();
		final int h = bi.getHeight();
		final Graphics2D g = bi.createGraphics();
      		double[] field = new double[2];
		g.clearRect(0, 0, w, h);
		for (int i = 0; i <= w; i += 21){
		for (int j = 0; j <= h; j += 21){ field = PDE.eField(c,i%21,j%21,z,l,W);
						  drawArrow(g, i, j, field[0]*10, field[1]);}}*/
		
	}

	// find max val and adjust for scale between 0,255 
	static int colourScale(double colourIn) {
		int colourOut = (int)(colourIn + 0.5d);
		return 5;
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
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {f.getGraphics().drawImage(bi, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}}, 0, 33);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {update(bi,W);}}}, 0, 3000);
	}

 	   static final double qPI = 0.25*Math.PI;
 	   static void drawArrow(final Graphics2D g, final int x1, final int y1, final double l1, final double a) {
		g.setColor(Color.WHITE);
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
