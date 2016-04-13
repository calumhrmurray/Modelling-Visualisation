import java.awt.Color;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;
import java.awt.Scrollbar;
import java.awt.*;


class PottsAnimate {	// -------------------------------------------------------------

  	static Label labelT = new Label("Temperature:	      ");

	static Ising Imodel = new Ising();

	static void init(final BufferedImage bi, final int W, final int c) {
			int[][] l = Imodel.initPotts(W,c);
			int val = 0;
			for (int i = 0; i < bi.getWidth(); i++){
			for (int j = 0; j < bi.getHeight(); j++){
			val = l[i][j];
			bi.setRGB(i, j, val == 0 ? Color.GREEN.getRGB() : 
					val == 1 ? Color.BLUE.getRGB() :
					val == 2 ? Color.RED.getRGB() : Color.WHITE.getRGB());}}
	}

	static void update(final BufferedImage bi, double t, final int d, final int W) {		
			// Work on l to update bi 
			int[][] l = Imodel.updatePotts(d,t,W);
			int val = 0;
			double M = Stats.updateM(l,W);
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){
			val = l[i][j];
			bi.setRGB(i, j, val == 0 ? Color.GREEN.getRGB() : 
					val == 1 ? Color.BLUE.getRGB() :
					val == 2 ? Color.RED.getRGB() : Color.WHITE.getRGB());}}
	}

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 5) throw new Exception("Arguments: width[pixels] height[pixels] period[milliseconds]");
		final int W = Integer.parseInt(args[0]);
		final int H = Integer.parseInt(args[0]);
		final double t = Double.parseDouble(args[1]); // initial temperature
		final int c = Integer.parseInt(args[2]);
		final int d = Integer.parseInt(args[3]);
		final int sim = Integer.parseInt(args[4]);

		final BufferedImage bg = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);	// background
		final BufferedImage fg = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);	// foreground
		final Frame f = new Frame();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate: Size -" + args[0] + " Init. Temp. - " + args[1] + " Initial -" + args[2] + " Dynamics -" + args[3]);
		f.setVisible(true);
		f.setSize(W, H + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});

		init(bg,W,c);
		fg.setData(bg.getData());

		new Timer().scheduleAtFixedRate(new TimerTask() {
			public void run() {f.getGraphics().drawImage(
				fg, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);
								}}, 0, 33);

		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {update(bg,t,d,W); fg.setData(bg.getData());}}, 0, 1);}
		//new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {init(bg,d,W); fg.setData(bg.getData());}}, 0, 1);}
}

	//----------------------------------------------------------------------




		
