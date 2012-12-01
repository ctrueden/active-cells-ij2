package ch.epfl.big.ij2.activecells;

import imagej.command.Command;
import imagej.data.Dataset;
import imagej.data.display.ImageDisplay;
import imagej.data.display.OverlayService;
import imagej.data.overlay.Overlay;
import imagej.data.overlay.PolygonOverlay;
import imagej.display.DisplayService;
import imagej.log.LogService;
import imagej.menu.MenuConstants;
import imagej.plugin.Menu;
import imagej.plugin.Parameter;
import imagej.plugin.Plugin;

import java.awt.Polygon;
import java.util.Collections;

import net.imglib2.RealPoint;
import net.imglib2.img.Img;
import net.imglib2.img.ImgPlus;
import net.imglib2.roi.PolygonRegionOfInterest;
import net.imglib2.roi.RegionOfInterest;
import net.imglib2.type.numeric.RealType;
import big.ij.snake2D.Snake2DKeeper;
import ch.epfl.big.ij2.activecells.snake.ESnake;
import ch.epfl.big.ij2.activecells.snake.ESnakeEnergyType;
import ch.epfl.big.ij2.activecells.snake.ESnakeParameters;
import ch.epfl.big.ij2.activecells.snake.ESnakeTargetType;
import ch.epfl.big.ij2.activecells.utils.ContourTracing;

@Plugin(label = "ESnake...", iconPath = "/icons/plugins/information.png", menu = {
		@Menu(label = MenuConstants.PLUGINS_LABEL),
		@Menu(label = "ESnake...", weight = 43) }, headless = true)
public class BIGSnake implements Command {

	@Parameter
	private Dataset dataSet;
	@Parameter
	private DisplayService displayService;
	@Parameter
	private OverlayService overlayService;
	@Parameter
	private ImageDisplay imageDisplay;
	@Parameter
	private LogService log;

	// ============================================================================
	// PUBLIC METHODS

	@Override
	public void run() {

		// Get input data
		ImgPlus<? extends RealType<?>> imgPlus = dataSet.getImgPlus();
		Img<? extends RealType<?>> inputImage = imgPlus.getImg();
		Overlay activeOverlay = overlayService.getActiveOverlay(imageDisplay);

		if (activeOverlay == null) {
			log.error("No ROI motherfucka!");
			return;
		}

		RegionOfInterest inputRegionOfInterest = activeOverlay
				.getRegionOfInterest();
		Polygon parsedInputRegionOfInterest = parseRegionOfInterest(inputRegionOfInterest);

		// Create snake objects

		// Snake parameters
		int maxLife = 500;
		int M = 3;
		double alpha = 0;
		boolean immortal = true;
		ESnakeTargetType detectType = ESnakeTargetType.BRIGHT;
		ESnakeEnergyType energyType = ESnakeEnergyType.REGION;
		ESnakeParameters eSnakeParameters = new ESnakeParameters(maxLife, M,
				alpha, immortal, detectType, energyType);

		// Image-related parameters
		double std = 5;

		ESnake mysnake = new ESnake(inputImage, std, eSnakeParameters,
				parsedInputRegionOfInterest, log);
		Snake2DKeeper keeper = new Snake2DKeeper();

		// Optimize the snake
		keeper.optimize(mysnake, null);
		// keeper.interactAndOptimize(mysnake, imp_);

		// Do something about the output
		setSnakeAsRegionOfInterest(mysnake);
	}

	// ============================================================================
	// PUBLIC METHODS

	private void setSnakeAsRegionOfInterest(ESnake snake) {
		// Removing the old ROI
		overlayService.removeOverlay(overlayService
				.getActiveOverlay(imageDisplay));

		PolygonOverlay newOverlay = new PolygonOverlay(dataSet.getContext());
		PolygonRegionOfInterest newRoi = newOverlay.getRegionOfInterest();
		for (int i = 0; i < snake.xPosSkin_.length; i++) {
			newRoi.addVertex(i, new RealPoint(snake.xPosSkin_[i],
					snake.yPosSkin_[i]));
		}
		// List<Overlay> overlays = overlayService.getOverlays(imageDisplay);
		// overlays.add(newOverlay);
		overlayService.addOverlays(imageDisplay,
				Collections.singletonList((Overlay) newOverlay));
		imageDisplay.update();
	}

	// ----------------------------------------------------------------------------

	private Polygon parseRegionOfInterest(RegionOfInterest roi) {

		boolean[] bitmask = new boolean[(int) ((roi.realMax(0) + 1) * ((roi
				.realMax(1) + 1)))];
		for (int x = (int) roi.realMin(0); x <= roi.realMax(0); x++) {
			for (int y = (int) roi.realMin(1); y <= (int) roi.realMax(1); y++) {
				bitmask[(int) (x + y * (roi.realMax(0) + 1))] = roi
						.contains(new double[] { x, y });
			}
		}

		ContourTracing tracer = new ContourTracing(bitmask,
				(int) (roi.realMax(0) + 1), (int) (roi.realMax(1) + 1));
		tracer.trace();
		return tracer.getTrace();
	}
}
