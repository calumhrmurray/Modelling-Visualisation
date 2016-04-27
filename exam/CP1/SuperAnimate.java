import java.awt.Color;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;

class SuperAnimate {	// -------------------------------------------------------------

	static Super Imodel = new Super();

	static void init(final BufferedImage bi, final int W, final double c) {
			int[][] l = Imodel.init(W,c);
			for (int i = 0; i < bi.getWidth(); i++){
			for (int j = 0; j < bi.getHeight(); j++){
			bi.setRGB(i, j, l[i][j]== -1 ? Color.RED.getRGB() : 
					l[i][j]== 0 ? Color.GREEN.getRGB() :
					l[i][j]== 1 ? Color.BLUE.getRGB() : Color.WHITE.getRGB());}}	
	}

	static void update(final BufferedImage bi, final int W, final double T) {
			int[][] l = Imodel.update(T,W);
			for (int i = 0; i < bi.getWidth(); i++){
			for (int j = 0; j < bi.getHeight(); j++){
			bi.setRGB(i, j, l[i][j]== -1 ? Color.RED.getRGB() : 
					l[i][j]== 0 ? Color.GREEN.getRGB() :
					l[i][j]== 1 ? Color.BLUE.getRGB() : Color.WHITE.getRGB());}}
	}
	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 4) throw new Exception("Arguments: width[pixels] height[pixels] T[]");
		final int W = Integer.parseInt(args[0]); // lattice size
		final double T = Double.parseDouble(args[1]); // initial temperature
		final double c = Double.parseDouble(args[2]); // initialisation
		final int sim = Integer.parseInt(args[3]); // simulation

		int H = W;

		final BufferedImage bi = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
		final Object lock = new Object();
		final Frame f = new Frame();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate " + args[0] + " " + args[1] + " " + args[2]);
		f.setVisible(true);
		f.setSize(W, H + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});

		init(bi,W,c);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {f.getGraphics().drawImage(bi, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}}, 0, 33);
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {synchronized(lock) {update(bi,W,T);}}}, 0, 33);
	}
}
