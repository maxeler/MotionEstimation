package motionEstimation;


import com.maxeler.maxcompiler.v2.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v2.kernelcompiler.KernelLib;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

/**
 * This class abstracts all the scalar inputs and other scalar constants based on the inputs that are used in other classes.
 * The scalars are accessible through public member fields.
 */
public class CoarseMeScalars extends KernelLib {

	public final static int windowX = 64;
	public final static int windowY = 32;

	public final static int maxPixelWidth  = 3840;
	public final static int maxPixelHeight = 2176;

	public final DFEVar widthInBlocks;
	public final DFEVar heightInBlocks;
	public final DFEVar widthPlusHaloInBlocks;
	public final DFEVar heightPlusHaloInBlocks;
	public final DFEVar numBlocks;
	public final DFEVar blocksToFill;
	public final DFEVar blocksToComplete;

	public final int haloBlocksX;
	public final int haloBlocksY;
	public final int blockSize;
	public final int numPipes;
	public final int maxWidthInBlocks;
	public final int maxHeightInBlocks;
	public final int maxBlocks;
	public final int maxBlocksToFill;
	public final int bufferDepth;


	/**
	 * This constructs an instance of scalar inputs. It does some Run Time checks on the inputs, it computes the values for member fields
	 * and creates the necessary scalar i/o's.
	 * @param owner      The Kernel that creates an instance of CourseMeScalars
	 * @param blockSize  the length of a side of the blocks in pixels
	 * @param numPipes   the number of pipes to be created in the DFE
	 */
	public CoarseMeScalars(Kernel owner, int blockSize, int numPipes) {
		super(owner);

		if (numPipes > blockSize) {
			throw new RuntimeException("Number of pipes can't be larger than block size.");
		}
		if (blockSize % numPipes != 0) {
			throw new RuntimeException("Number of pipes must be a factor of block size.");
		}
		if (!MathUtils.isPowerOf2(blockSize)) {
			throw new RuntimeException("Block size must be a power of 2.");
		}
		if (numPipes <=0 || blockSize <= 0) {
			throw new RuntimeException("Number of pipes and block size must both be positive integers.");
		}

		haloBlocksX = MathUtils.ceilDivide(windowX/2, blockSize);
		haloBlocksY = MathUtils.ceilDivide(windowY/2, blockSize);
		maxWidthInBlocks  = (maxPixelWidth / 2 + windowX) / blockSize;
		maxHeightInBlocks = (maxPixelHeight / 2) / blockSize;
		maxBlocks         = maxWidthInBlocks * maxHeightInBlocks;
		maxBlocksToFill   = maxWidthInBlocks * (2 * haloBlocksY + 1);
		bufferDepth       = MathUtils.nextPowerOfTwo(maxBlocksToFill * blockSize / numPipes);

		this.blockSize = blockSize;
		this.numPipes  = numPipes;

		widthInBlocks          = io.scalarInput("widthInBlocks",         dfeUInt(MathUtils.bitsToRepresent(maxWidthInBlocks)));
		heightInBlocks         = io.scalarInput("heightInBlocks",        dfeUInt(MathUtils.bitsToRepresent(maxHeightInBlocks)));
		widthPlusHaloInBlocks  = io.scalarInput("widthPlusHaloInBlocks", dfeUInt(MathUtils.bitsToRepresent(maxWidthInBlocks)));
		heightPlusHaloInBlocks  = io.scalarInput("heightPlusHaloInBlocks", dfeUInt(MathUtils.bitsToRepresent(maxHeightInBlocks)));

		DFEType blockCountType = dfeUInt(MathUtils.bitsToRepresentUnsigned(maxBlocks));
		numBlocks        = io.scalarInput("numBlocks",        blockCountType);
		blocksToFill     = io.scalarInput("fillBlocks",       blockCountType);
		blocksToComplete = io.scalarInput("blocksToComplete", blockCountType);

	}
}