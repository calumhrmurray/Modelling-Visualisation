import java.awt.Color;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;
import java.io.*;
import java.awt.*;
import java.awt.event.*;



class Animate2 {	// -------------------------------------------------------------

	static boolean running = true;               
	static Sirs Smodel = new Sirs();

	static int green = Color.MAGENTA.getRGB();
	static int black = Color.PINK.getRGB();
	static int grey = Color.WHITE.getRGB();
	static int white = Color.WHITE.getRGB();

	static void init(final BufferedImage bi, final int W) {
			int[][] l = Smodel.init(W);
			for (int i = 0; i < bi.getWidth(); i++){
			for (int j = 0; j < bi.getHeight(); j++){
			bi.setRGB(i, j, l[i][j]== -1 ? Color.GREEN.getRGB() : 
					l[i][j]== 0 ? Color.GRAY.getRGB() :
					l[i][j]== 1 ? Color.DARK_GRAY.getRGB() : Color.WHITE.getRGB());}}

	}


	
	static void update(final BufferedImage bi, final int W,final double p1, final double p2, final double p3) {		
		int[][] l = Smodel.update(W,p1,p2,p3); 
			
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){
			int a = l[i][j]; 
			bi.setRGB(i, j, a == -1 ? green : 
					a == 0 ? black :
					a == 1 ? grey : white);}}
	}


	

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 4) throw new Exception("Arguments: width[pixels] height[pixels] period[milliseconds]");
		final int W = Integer.parseInt(args[0]);
		final double p1 = Double.parseDouble(args[1]);
		final double p2 = Double.parseDouble(args[2]);
		final double p3 = Double.parseDouble(args[3]);

		final BufferedImage bi = new BufferedImage(W, W, BufferedImage.TYPE_INT_RGB);
		final Object lock = new Object();
		final Frame f = new Frame();
		Panel controlPanel = new Panel();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate " + args[0]);
		f.setVisible(true);
		f.setSize(W, W + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});


		init(bi,W);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {f.getGraphics().drawImage(bi, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}}, 0, 33);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {
				if (running) {update(bi,W,p1,p2,p3);}}}}, 0, 100);

	}

	//----------------------------------------------------------------------

}


