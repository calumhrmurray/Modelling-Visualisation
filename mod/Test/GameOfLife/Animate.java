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


class Animate {	// -------------------------------------------------------------

  	static Label labelT = new Label("Temperature:	      ");
	
	static Game Gmodel = new Game();

	static void init(final BufferedImage bi, final int W, final int c) {
			int[][] l = Gmodel.init(W,c);
			for (int i = 0; i < bi.getWidth(); i++){
			for (int j = 0; j < bi.getHeight(); j++){
			bi.setRGB(i, j, l[i][j]==1 ? Color.GREEN.getRGB() : Color.WHITE.getRGB());}}	
	}

	static void update(final BufferedImage bi, final int W) {		
			int[][] l = Gmodel.update(W);
			for (int i = 0; i < W; i++){
			for (int j = 0; j < W; j++){
			bi.setRGB(i, j, l[i][j]==1 ? Color.GREEN.getRGB() : Color.WHITE.getRGB());}}		
	}

	//----------------------------------------------------------------------

	public static void main(final String[] args) throws Exception {
		if (args.length != 2) throw new Exception("Arguments: width[pixels] height[pixels] period[milliseconds]");
		final int W = Integer.parseInt(args[0]);
		final int c = Integer.parseInt(args[1]);


		final BufferedImage bg = new BufferedImage(W, W, BufferedImage.TYPE_INT_RGB);	// background
		final BufferedImage fg = new BufferedImage(W, W, BufferedImage.TYPE_INT_RGB);	// foreground
		final Frame f = new Frame();
		f.setIgnoreRepaint(true);
		f.setTitle("Animate: Size -" + args[0]);
		f.setVisible(true);
		f.setSize(W, W + f.getInsets().top);
		f.setExtendedState(Frame.MAXIMIZED_BOTH);
		f.addWindowListener(new WindowAdapter() {public void windowClosing(WindowEvent we) {System.exit(0);}});

			
		init(bg,W,c);
		fg.setData(bg.getData());

		new Timer().scheduleAtFixedRate(new TimerTask() {
			public void run() {f.getGraphics().drawImage(fg, 0, f.getInsets().top, f.getWidth(), f.getHeight() - f.getInsets().top, null);}}, 0, 1000);
	
		new Timer().scheduleAtFixedRate(new TimerTask() {public void run() {update(bg,W); fg.setData(bg.getData());}}, 0, 1000);}

}

	//----------------------------------------------------------------------




		
