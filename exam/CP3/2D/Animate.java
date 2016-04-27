import java.awt.Color;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;

class Animate {	// -------------------------------------------------------------

	static PDE Pmodel = new PDE();

	static void init(final BufferedImage bi, final int W) {
		double[][] c = Pmodel.init(W);
		int colour = 0;
		// print image	
		for (int x = 0; x < bi.getWidth(); x++){
		for (int y = 0; y < bi.getHeight(); y++){ colour = (int)(c[x][y]+0.5d); 
							  bi.setRGB(x, y, new Color(colour,colour, colour ).getRGB());}}
	}

	static void update(final BufferedImage bi, final int W) {
		double[][] c = Pmodel.update(W); 
		int colour = 0;
		// print image	
		for (int x = 0; x < bi.getWidth(); x++){
		for (int y = 0; y < bi.getHeight(); y++){ colour = (int)(c[x][y]+0.5d); 
							  bi.setRGB(x, y, new Color(colour,colour, colour ).getRGB());}}
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
		

		final BufferedImage bi = new BufferedImage(W, W, BufferedImage.TYPE_INT_RGB);
		final Object lock = new Object();
		final Frame f = new Frame();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate " + args[0] + " " + args[1] + " " + args[2]);
		f.setVisible(true);
		f.setSize(W, W + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});

		init(bi,W);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {f.getGraphics().drawImage(bi, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}}, 0, 33);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {update(bi,W);}}}, 0, 20);
	}
}
