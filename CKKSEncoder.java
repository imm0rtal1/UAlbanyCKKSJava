/* https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/ was used for any Matrix/Complex/Polynomials */

import java.util.Random;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
/*https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/linear/RealMatrix.html */
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.LUDecomposition;
public class CKKSEncoder {
	/**
	 * Power of 2 for the matrix, with half of M being the size of stored data
	 */
	private int M;
	/**
	 * secret key for the matrix
	 */
	private Complex xi;
	/**
	 * Public key for the matrix, of size M/2 by M/2
	 */
	private RealMatrix basis;
	/**
	 * scale variable for the matrix
	 */

	private float scale;
	/**
	 * Initializes a CKKS encoder
	 * @param m M value for the encoder
	 * @param scaleVal Scale value for the encoder
	 */
	public CKKSEncoder(int m, float scaleVal) {
		M=m;
		xi= Complex.valueOf(Math.exp((Complex.I).add(1).multiply(2*Math.PI/m).getReal()));
		basis = Create_sigma_R_Basis(xi,m);
		scale=scaleVal;
	}
	@SuppressWarnings("null")
	/**
	 * creates the vandermonde matrix for the encoder 
	 * @param xi the xi value for the matrix
	 * @param M the M value for the matrix
	 * @return the vandermonde matrix
	 * @throws NotStrictlyPositiveException
	 */
	public static RealMatrix Vandermonde(Complex xi, int M) throws NotStrictlyPositiveException {
		int N = M/2;
		RealMatrix matrix = null;
		matrix = matrix.createMatrix(N, N);
		for(int i = 0; i < N; i++){
			Complex root = xi.pow(2*(i+1));
			RealVector row = null;
			for(int j = 0; j<N;j++) {
				Complex value =root.pow(j);
				row.append(value.getReal());
			}
			matrix.setRowVector(i,row);
		}
		return matrix;
	}
	/**
	 * encodes vector b into a polynomial
	 * @param b the vector being encoded
	 * @return the polynomial ciphertext produced
	 */
	public PolynomialFunction sigmaInverse(RealVector b) {
		RealMatrix A = CKKSEncoder.Vandermonde(xi,M);
		LUDecomposition decomposition = new LUDecomposition(A);
		RealVector coeff = decomposition.getSolver().solve(b);
		PolynomialFunction p = new PolynomialFunction(coeff.toArray());
		return p;
	}
	/**
	 * Decodes polynomial p into a vector
	 * @param p the polynomial being decoded
	 * @return the decoded vector produced
	 */
	public RealVector Sigma(PolynomialFunction p) {
		int N=M/2;
		RealVector outputs = null;
		outputs = outputs.map(p);
		for(int i=0;i<N;i++) {
			Complex root = xi.pow(2*(i+1));
			double output= p.value(root.getReal());
			outputs.setEntry(i, output);
		}
		return outputs;
	}
	/**
	 * determines difference between reconstructed vector and original one
	 * @param b original vector that will be encoded and then decoded
	 * @return error from encoding
	 */
	public double FindError(RealVector b) {
		PolynomialFunction p = this.sigmaInverse(b);
		RealVector B_reconstructed = this.Sigma(p);
		b.subtract(B_reconstructed);
		double norm = b.getNorm();
		return norm;	
	}
	/**
	 * adds two encoded values together
	 * @param p1 encoded value 1
	 * @param p2 encoded value 2
	 * @return sum of encoded values
	 */
	public static PolynomialFunction addPolynomials(PolynomialFunction p1,  PolynomialFunction p2) {
		return(p1.add(p2));
	}
	/**
	 * multiplies two encoded values together
	 * @param p1 encoded value 1
	 * @param p2 encoded value 2
	 * @return product of p1 and p2
	 */
	public static PolynomialFunction multiplyPolynomials(PolynomialFunction p1,  PolynomialFunction p2) {
		int N=p2.degree();
		double[] modArray= new double[N+1];
		modArray[N]=1;
		modArray[0]=1;
		for(int i=1;i<N;i++) {
			modArray[i]=0;
		}
		PolynomialFunction poly_modulo = new PolynomialFunction(modArray);
		PolynomialFunction p_mult = p1.multiply(p2);
		RealVector multV = null;
		multV = multV.map(p_mult);
		RealVector modV = null;
		modV = modV.map(poly_modulo);
		double[] tMod = modV.toArray();
		double[] tMult = multV.toArray();
		double[] results = new double[N];
		for(int i=0;i<tMod.length;i++) {
			results[i]=tMult[i]% tMod[i];
		}
		p_mult = new PolynomialFunction(results);
		return p_mult;
	}
	/**
	 * Converts array of size M/2 to one of size M/4
	 * @param z array of size M/2
	 * @return array of size M/4
	 */
	public Complex[] pi(Complex[] z) {
		int N=M/4;
		Complex[] output = new Complex[N];
		for(int i=0;i<N;i++) {
			output[i]=z[i];
		}
		return output;
	}
	/**
	 * converts array of size M/4 to one of size M/2
	 * @param z array of size M/4
	 * @return array of size M/2
	 */
	public Complex[] pi_inverse(Complex[] z) {
		int N=M/2;
		Complex[] output= new Complex[N];
		for(int i=0;i<N/2;i++) {
			output[i]=z[i];
		}
		Complex[] z_conjugate = new Complex[N/2];
		for(int i=0;i<N/2;i++) {
			z_conjugate[i]=z[z.length-1-i];
			z_conjugate[i]=z_conjugate[i].conjugate();
		}
		for(int i=0;i<N/2;i++) {
			int placement = i+N;
			output[placement]=z_conjugate[i];
		}
		return output;
	}
	/**
	 * creates basis matrix for the encoder
	 * @param xi the xi value for the basis
	 * @param M M value for the basis
	 * @return the basis matrix for the encoder
	 */
	public static RealMatrix Create_sigma_R_Basis(Complex xi, int M) {
		RealMatrix output = Vandermonde(xi,M);
		return(output.transpose());
	}
	@SuppressWarnings("null")
	/**
	 * calculates coords for computing basis coords 
	 * @param coords default coords not yet in basis
	 * @return coords in basis
	 */
	public double[] calculateB(double[] coords) {
		RealMatrix basis_t=basis.transpose();
		RealMatrix coordsMatrix = null;
		RealVector v = null;
		for(int i=0;i<coords.length;i++) {
			v.append(coords[i]);
		}
		coordsMatrix=coordsMatrix.createMatrix(coords.length, 1);
		coordsMatrix.setColumnVector(0, v);
		return basis_t.multiply(coordsMatrix).getColumn(0);
	}
	/**
	 * computes the coords from the matrix 
	 * @param Z the array of coords being translated to the matrix
	 * @return the translated coords
	 */
	public double[] compute_basis_coordinates(Complex[] Z) {
		double[] returnArray = new double[M/2];
		//Convert Z into a RealVector
		RealVector v = null;
		for(int i=0;i<Z.length;i++) {
			v.setEntry(i,Z[i].getReal());
		}
		//calculate basis coords
		for(int i=0;i<Z.length;i++) {
			 RealVector b = basis.getColumnVector(i);
			 double VdotB = v.dotProduct(b);
			 double VdotV = v.dotProduct(v);
			 returnArray[i]=VdotB/VdotV;
		}
		return returnArray;
	}
	/**
	 * rounds a set of coordinates
	 * @param coords the coordinates being rounded
	 * @return the rounded coordinates
	 */
	public static double[] round_coordinates(double[] coords) {
		for (int i=0;i<coords.length;i++) {
			coords[i]=coords[i]-Math.floor(coords[i]);
		}
		return coords;
	}
	/**
	 * randomly rounds coordinates 
	 * @param coords the coordinates being randomly rounded 
	 * @return the randomly rounded coordinates
	 */
	public static int[] coordinate_wise_random_rounding(double[] coords) {
		double[] r = round_coordinates(coords);
		Random random = new Random();
		int[] f = new int[r.length];
		for(int i=0;i<r.length;i++) {
			double rVal = 1-r[i];
			double rand = random.nextDouble();
			if(rVal>rand) {
				f[i]=(int)(coords[i]-1);
			}
			else {
				f[i]=(int)coords[i];
			}
		}
		int[] rounded_coords = new int[coords.length];
		for (int i=0;i<rounded_coords.length;i++) {
			rounded_coords[i]=(int)(coords[i]-f[i]);
		}
		return rounded_coords;
	}
	/**
	 * converts coordinates into a discrete set of coordinates for encoding
	 * @param z the array of coordinates being translated 
	 * @return the discrete set of coordinates used for encoding
	 */
	@SuppressWarnings("null")
	public double[] sigma_R_discretization(Complex[] z) {
		double[] coords= compute_basis_coordinates(z);
		int[] rounded_coords = coordinate_wise_random_rounding(coords);
		RealMatrix m = null;
		m.createMatrix(z.length, 1);
		for(int i=0;i<z.length;i++) {
			coords[i]=(double)rounded_coords[i];
		}
		m.setColumn(0, coords);
		RealMatrix y = basis.transpose().multiply(m);
		return y.getColumn(0);
	}
	
	/**
	 * Encodes an array of complex numbers into a polynomial function
	 * @param z the array of complex numbers being encoded
	 * @return the encoded polynomial function
	 */
	@SuppressWarnings("null")
	public PolynomialFunction encode(Complex[] z) {
		Complex[] pi_z = this.pi_inverse(z);
		for(int i=0;i<pi_z.length;i++) {
			pi_z[i]=pi_z[i].multiply(scale);
		}
		double[] rounded_pi_z = sigma_R_discretization(pi_z);
		//convert rounded_pi_z to realVector roundedPiZ
		RealVector roundedPiZ = null;
		for(int i=0;i<pi_z.length;i++) {
			roundedPiZ.setEntry(i, rounded_pi_z[i]);
		}
		PolynomialFunction p = sigmaInverse(roundedPiZ);
		int[] coef = new int[pi_z.length];
		double[] coeffs = p.getCoefficients();
		for(int i=0;i<pi_z.length;i++) {
			coef[i]=(int)coeffs[i];
			coeffs[i]=coef[i];
		}
		p=new PolynomialFunction(coeffs);
		return p;
	}
	/**
	 * decodes polynomial function P into an array of complex numbers
	 * @param p the polynomial function being decoded
	 * @return the decoded array of complex numbers
	 */
	public Complex[] decode(PolynomialFunction p) {
		double[] pCoeffs = p.getCoefficients();
		for (int i=0;i<pCoeffs.length;i++) {
			pCoeffs[i]=pCoeffs[i]/scale;
		}
		p = new PolynomialFunction(pCoeffs);
		double[] zDouble = Sigma(p).toArray();
		Complex[] z = new Complex[pCoeffs.length];
		for(int i=0;i<pCoeffs.length;i++) {
			z[i]=new Complex(zDouble[i]);
		}
		Complex[] Pi_z = pi(z);
		return Pi_z;
	}
}
//To Test: take array of size M/2, use calculateB function on it, then use that value for encode/decode functions