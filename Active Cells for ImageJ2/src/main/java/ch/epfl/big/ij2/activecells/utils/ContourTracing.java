/*
Copyright (c) 2010-2012 Thomas Schaffter & Ricard Delgado-Gonzalo

WingJ is licensed under a
Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported License.

You should have received a copy of the license along with this
work. If not, see http://creativecommons.org/licenses/by-nc-nd/3.0/.

If this software was useful for your scientific work, please cite our paper(s)
listed on http://lis.epfl.ch/wingj.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

package ch.epfl.big.ij2.activecells.utils;

import ij.IJ;

import java.awt.Point;
import java.awt.Polygon;
import java.util.Vector;

/** 
 * Square tracing algorithm.
 * <p>
 * The idea behind the square tracing algorithm is very simple; this could be 
 * attributed to the fact that the algorithm was one of the first attempts to 
 * extract the contour of a binary pattern.
 * <p>
 * Source: http://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/square.html
 * 
 * @version August 12, 2011
 * 
 * @author Ricard Delgado-Gonzalo (ricard.delg...@gmail.com)
 */
public class ContourTracing {
	
	private boolean[] T;
	private int width;
	private int height;
	private Vector<Point> path;
	private static final int UP = 0;
	private static final int DOWN = 1;
	private static final int LEFT = 2;
	private static final int RIGHT = 3;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor. */
	public ContourTracing(boolean[] T, int width, int height){
		this.T = T;
		this.width = width;
		this.height = height;
		path = new Vector<Point>();
	}
	
	// ----------------------------------------------------------------------------

	/** Constructor. */
	public ContourTracing(float[] T, int width, int height){
		this.T = new boolean[width*height];
		
   		for(int i =0; i<T.length; i++){
   			if(T[i]!=0){
   				this.T[i] = true;
   			}else{
   				this.T[i] = false;
   			}
   		}
		
		this.width = width;
		this.height = height;
		path = new Vector<Point>();
	}
	
	// ----------------------------------------------------------------------------

	/** Traces the outline of a binary region using the Square Tracing Algorithm. */
	public void trace() {
		int sx = 0;
		int sy = 0;
		int px = 0;
		int py = 0;
		int direction = UP;
		
		boolean found = false;
		for(int y=height-1; y>=0 && !found; y--){
			for(int x=0; x<width && !found; x++){
				if(T[x+y*width]){
					found = true;
					sx = x;
					sy = y;
					px = x;
					py = y;
					path.addElement(new Point(x,y));
				}
			}
		}
		
		direction = RIGHT;
		px--;
		boolean cond;
		while(!(sx==px && sy==py)){
			
			if(px<0 || px>=width || py<0 || py>=height){
				cond=false;
			}else{
				cond = T[px+py*width]; 
			}
			if(cond){
				path.addElement(new Point(px,py));
				switch(direction){
					case UP:
						px++;
						direction = LEFT;
						break;
					case DOWN:
						px--;
						direction = RIGHT;
						break;
					case LEFT:
						py--;
						direction = DOWN;
						break;
					case RIGHT: 
						py++;
						direction = UP;
						break;
					default:
						IJ.error("Orientation unknown.");
						break;
				}
			}else{
				switch(direction){
				case UP:
					px--;
					direction = RIGHT;
					break;
				case DOWN:
					px++;
					direction = LEFT;
					break;
				case LEFT:
					py++;
					direction = UP;
					break;
				case RIGHT: 
					py--;
					direction = DOWN;
					break;
				default:
					IJ.error("Orientation unknown.");
					break;
				}
			}
		}
		eliminateDuplicates();
	}

	// ============================================================================
	// PRIVATE METHODS
	
	/** Removes duplicated points. */
	private void eliminateDuplicates() {
		if(path!=null){
			if(path.size()>1){
				for(int i=0; i<path.size()-1; i++){
					if(path.elementAt(i).x == path.elementAt(i+1).x && path.elementAt(i).y == path.elementAt(i+1).y){
						path.removeElementAt(i+1);
						i--;
					}
				}
			}
		}
	}

	// ============================================================================
	// SETTERS AND GETTERS

	/** Retrieves the trace as a Polygon. */
	public Polygon getTrace(){
		return(new Polygon(getXCoordinates(), getYCoordinates(), getNPoints()));
	}

	// ----------------------------------------------------------------------------

	/** Retrieves the x's coordinates of the trace. */
	public int[] getXCoordinates() {
		int nPoints = path.size();
		int[] xCoordinates = new int[nPoints];
		for(int i=0; i<path.size(); i++){
			xCoordinates[i] = path.elementAt(i).x;
		}
		return xCoordinates;
	}
	
	// ----------------------------------------------------------------------------

	/** Retrieves the y's coordinates of the trace. */
	public int[] getYCoordinates() {
		int nPoints = path.size();
		int[] yCoordinates = new int[nPoints];
		for(int i=0; i<nPoints; i++){
			yCoordinates[i] = path.elementAt(i).y;
		}
		return yCoordinates;
	}	
	
	// ----------------------------------------------------------------------------
	
	/** Retrieves the total number of points in the trace. */
	public int getNPoints() {
		return 	path.size();
	}
}