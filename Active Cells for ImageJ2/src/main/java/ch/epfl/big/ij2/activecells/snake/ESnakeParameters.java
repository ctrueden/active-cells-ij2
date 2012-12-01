package ch.epfl.big.ij2.activecells.snake;

/**
 * Class that wraps the parameters of the E-Snake.
 * 
 * @version November 13, 2012
 * 
 * @author Ricard Delgado-Gonzalo (ricard.delgado@gmail.com)
 */
public class ESnakeParameters {

	/**
	 * If <code>true</code> indicates that the snake will keep iterating till
	 * the optimizer decides so.
	 */
	private boolean immortal_ = true;

	/** Number of spline vector coefficients. */
	private int M_ = 0;

	/** Energy tradeoff factor. */
	private double alpha_ = 0;

	/** Maximum number of iterations allowed when the immortal is false. */
	private int maxLife_ = 0;

	/** Indicates the energy function of the snake. **/
	private ESnakeEnergyType energyType_;

	/** Indicates the type of features to detect (bright or dark), **/
	private ESnakeTargetType detectType_;

	// ============================================================================
	// PUBLIC METHODS

	/** Default constructor. */
	public ESnakeParameters(int maxLife, int M, double alpha, boolean immortal,
			ESnakeTargetType detectType, ESnakeEnergyType energyType) {
		setMaxLife(maxLife);
		setM(M);
		setAlpha(alpha);
		setImmortal(immortal);
		setDetectType(detectType);
		setEnergyType(energyType);
	}

	// ----------------------------------------------------------------------------

	@Override
	public String toString() {
		return new String("[E-Snake parameters: immortal = " + immortal_
				+ ", M = " + M_ + ", alpha = " + alpha_
				+ ", maximum number of iterations = " + maxLife_
				+ ", energy type = " + energyType_ + ", detect type = "
				+ detectType_ + "]");
	}

	// ============================================================================
	// GETTERS AND SETTERS

	public double getAlpha() {
		return alpha_;
	}

	// ----------------------------------------------------------------------------

	public void setAlpha(double alpha) {
		alpha_ = alpha;
	}

	// ----------------------------------------------------------------------------

	public boolean isImmortal() {
		return immortal_;
	}

	// ----------------------------------------------------------------------------

	public void setImmortal(boolean immortal) {
		immortal_ = immortal;
	}

	// ----------------------------------------------------------------------------

	public int getM() {
		return M_;
	}

	// ----------------------------------------------------------------------------

	public void setM(int M) {
		M_ = M;
	}

	// ----------------------------------------------------------------------------

	public ESnakeTargetType getDetectType() {
		return detectType_;
	}

	// ----------------------------------------------------------------------------

	public void setDetectType(ESnakeTargetType detectType) {
		detectType_ = detectType;
	}

	// ----------------------------------------------------------------------------

	public ESnakeEnergyType getEnergyType() {
		return energyType_;
	}

	// ----------------------------------------------------------------------------

	public void setEnergyType(ESnakeEnergyType energyType) {
		energyType_ = energyType;
	}

	// ----------------------------------------------------------------------------

	public int getMaxLife() {
		return maxLife_;
	}

	// ----------------------------------------------------------------------------

	public void setMaxLife(int maxLife) {
		maxLife_ = maxLife;
	}
}
