package motionEstimation;


import static com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEStructType.sft;

import java.util.ArrayList;
import java.util.List;

import com.maxeler.maxcompiler.v2.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v2.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.Reductions;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.CounterChain;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEStruct;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEStructType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEStructType.StructFieldType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVector;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVectorType;
import com.maxeler.maxcompiler.v2.utils.Bits;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

/**
 * This kernel performs the calculations to do a full block matching motion estimation.
 * It read the source and reference frames as inputs and provides the motion vectors as outputs.
 * The frame size, the block size, and number of pipes could be changed through parameters or
 * scalar inputs.
 * @author skamali
 *
 */
public class CoarseMeKernel extends Kernel {

	private final CoarseMeScalars m_scalars;

	/**
	 * Constructs the kernel
	 * @param parameters  Kernel parameters
	 * @param blockSize   the length of a block side in pixels
	 * @param numPipes    number of pipes to be constructed. More pipes provides better performance in expense of computation area
	 * @param inputBits   The width of input elements in bits
	 */
	protected CoarseMeKernel(KernelParameters parameters, int blockSize, int numPipes, int inputBits) {
		super(parameters);
		m_scalars = new CoarseMeScalars(this, blockSize, numPipes);
		CoarseMeBuffer refBuffer = new CoarseMeBuffer(m_scalars);  // creates the buffer to hold the reference frame blocks
		DFEVar isRunning = refBuffer.getOutputEnable();

		CounterChain chain = control.count.makeCounterChain(isRunning);
		chain.addCounter(m_scalars.heightInBlocks, 1);
		chain.addCounter(m_scalars.widthInBlocks, 1);
		DFEVar row = chain.addCounter(CoarseMeScalars.windowY, 1);
		DFEVar col = chain.addCounter(CoarseMeScalars.windowX, numPipes);

		DFEVar startOfBlock = (row === 0) & (col === 0) & isRunning;
		DFEVar endOfBlock   = (row === CoarseMeScalars.windowY - 1) & (col === CoarseMeScalars.windowX - numPipes);
		DFEVar[][] source = getSource(startOfBlock, blockSize, numPipes, inputBits);


		List<DFEVar[][]> sads = new ArrayList<DFEVar[][]>(numPipes);
		for (int pipe = 0; pipe < numPipes; pipe++) {
			DFEVar[][] reference = refBuffer.getReferenceBlock(pipe);
			sads.add(computeSad(source, reference));

		}

		DFEVar[][] vectorsAndCosts = calcMotionVectors(row, col);

		DFEVar[][][] minimumThisCycle = getCycleMinimum(sads, vectorsAndCosts);

		DFEVar[][][] minimum = getAggregateMinimum(startOfBlock, minimumThisCycle);


		DFEStruct output = getMotionVectors(minimum);
		io.output("mv", output.getType(), endOfBlock) <== output;

	}

	/**
	 * This function gets a source block from the input. It keeps the block in the pipeline until processing the current
	 * search window is completed. Then, the next source block is brought into the pipeline.
	 * @param startOfBlock  a signal that signifies that processing a new source block is just started
	 * @param blockSize     the length of a block side in pixels
	 * @param numPipes      the number of pipes
	 * @param inputBits     the input width in bits
	 * @return
	 */
	private DFEVar[][] getSource(DFEVar startOfBlock, int blockSize, int numPipes, int inputBits) {
		DFEVar[][] output = new DFEVar[blockSize][blockSize];
		int cyclesToReadBlock = blockSize * blockSize / numPipes;
		DFEVar inputEnable = startOfBlock;

		// looks ahead to see when the source block is needed, as reading a source blocks takes more than one cycle
		for (int i = 1; i < cyclesToReadBlock; i++) {
			inputEnable = inputEnable # stream.offset(startOfBlock, i);
		}
		inputEnable = inputEnable !== 0;

		DFEVectorType<DFEVar> type = new DFEVectorType<DFEVar>(dfeUInt(10), numPipes);
		DFEVector<DFEVar> source = io.input("source", new DFEVectorType<DFEVar>(dfeUInt(inputBits), numPipes), inputEnable).cast(type);

		for (int i = 0; i < cyclesToReadBlock; i++) {
			DFEVector<DFEVar> temp = stream.offset(source, i + 1 - cyclesToReadBlock);
			for (int j = 0; j < numPipes; j++) {
				int row = (j + i * numPipes) / blockSize;
				int col = (j + i * numPipes) % blockSize;
				// holds the block in the pipeline until the search window processing is completed
				output[row][col] = Reductions.streamHold(temp[j], startOfBlock);
			}
		}

		return output;
	}

	/**
	 * @return  the absolute difference between the two inputs
	 */
	private DFEVar absDiff(DFEVar a, DFEVar b){
		optimization.pushEnableBitGrowth(true);
		DFEVar diff = a - b;
		optimization.popEnableBitGrowth();
		return diff < 0  ? -diff : diff;
	}

	/**
	 * a recursive function that produces an adder tree to add a number of different variables together
	 * @param input  the variables that need to be added together
	 * @return       the result of the add
	 */
	private DFEVar adderTree(DFEVar... input) {
		if(input.length == 0)
			throw new RuntimeException("adderTree called with empty list");
		if (input.length == 1) {
			return input[0];
		}

		DFEVar[] output = new DFEVar[MathUtils.ceilDivide(input.length, 2)];
		for (int i = 0; i < input.length / 2; i++) {
			output[i] = input[2 * i] + input[2 * i + 1];
		}
		if (input.length % 2 == 1) {
			output[input.length / 2] = input[input.length - 1];
		}

		return adderTree(output);
	}

	/**
	 * Computes the Sum of Absolute difference between the source block and the reference block, it computes the SAD for all the sub-blocks
	 * @param source     2D array holding the source block
	 * @param reference  2D array holding the reference block
	 * @return           The SADS organised as a 2D array as output[size of sub-block][sub-block index]
	 */
	private DFEVar[][] computeSad(DFEVar[][] source, DFEVar[][] reference) {
		int blockSize = source.length;
		DFEVar[][] output = new DFEVar[MathUtils.bitsToRepresent(blockSize / 4)][blockSize / 4 * blockSize / 4];
		DFEVar[] sadList = new DFEVar[16];

		// computing SAD for 4 x 4 sub-blocks
		for (int blockRow = 0; blockRow < blockSize / 4; blockRow++) {
			for (int blockCol = 0; blockCol < blockSize / 4; blockCol++) {
				List<DFEVar> srcBlock = new ArrayList<DFEVar>();
				List<DFEVar> refBlock = new ArrayList<DFEVar>();
				for (int row = 0; row < 4; row++) {
					for (int col = 0; col < 4; col++) {
						srcBlock.add(source[row + blockRow * 4][col + blockCol * 4]);
						refBlock.add(reference[row + blockRow * 4][col + blockCol * 4]);
						sadList[row*4 + col] = absDiff(source[row + blockRow * 4][col + blockCol * 4], reference[row + blockRow * 4][col + blockCol * 4]);
					}
				}
				output[0][blockCol + blockRow * blockSize / 4] = adderTree(sadList);
			}
		}

		// computing SAD for the rest of the sub-blocks
		int index = 1;
		for (int i = 8; i <= blockSize; i *= 2) {
			int numBlocks = (blockSize * blockSize) / (i * i);
			for (int j = 0; j < numBlocks; j++) {
				output[index][j] = sumSads(output[index - 1], j, blockSize / i);
			}
			index++;
		}

		return output;
	}

	/**
	 * This method computes SAD of a block by adding the SAD of the smaller sub-blocks within it
	 * @param sads           sads 2D array as organised by method comuteSAd()
	 * @param block          the index that iterates through all the sub-blocks of the same size
	 * @param widthInBlocks  width of the regular block in terms of the current sub-block size
	 * @return               outputs the sum
	 */
	private DFEVar sumSads(DFEVar[] sads, int block, int widthInBlocks) {
		int row = 2 * (block / widthInBlocks);
		int col = 2 * (block % widthInBlocks);
		optimization.pushEnableBitGrowth(true);
		DFEVar top    = sads[row       * 2 * widthInBlocks + col] + sads[row       * 2 * widthInBlocks + col + 1];
		DFEVar bottom = sads[(row + 1) * 2 * widthInBlocks + col] + sads[(row + 1) * 2 * widthInBlocks + col + 1];
		DFEVar output = top + bottom;
		optimization.popEnableBitGrowth();

		return output;
	}

	/**
	 *Returns the motion vector for the block being processed
	 * @param row  The row number with reference to the current active search window
	 * @param col  The column number with reference to the current active search window
	 * @return     returns the motion vectors as a 2d array, first dimension determines the pipe number, second dimension determines x or y
	 */
	private DFEVar[][] calcMotionVectors(DFEVar row, DFEVar col) {
		DFEVar[][] output = new DFEVar[m_scalars.numPipes][2];

		DFEVar currY = getCurrentY(row);
		DFEVar[] currX = getCurrentX(col, m_scalars.numPipes);

		optimization.pushEnableBitGrowth(true);

		for (int i = 0; i < m_scalars.numPipes; i++) {
			output[i][0] = currX[i];
			output[i][1] = currY;
		}
		optimization.popEnableBitGrowth();

		return output;
	}

	/**
	 * Gets the motion vector to the current reference blocks in the x direction
	 * @param col       the current column number within the search window
	 * @param numPipes  total number of pipes
	 * @return          output is in the form of 1D array for the number of pipes
	 */
	private DFEVar[] getCurrentX(DFEVar col, int numPipes) {
		DFEVar[] output = new DFEVar [numPipes];
		output[0] = 2 * col.cast(dfeInt(col.getType().getTotalBits() + 2)) - CoarseMeScalars.windowX;
		for (int i = 1; i < numPipes; i++) {
			output[i] = output[0] + 2 * i;  // the blocks for all the pipes are next to each other. A simple loop generates the X for the other pipes
		}

		return output;
	}

	/**
	 * Gets the motion vector to the current reference blocks in the y direction
	 * @param row   the current row number within the search window
	 */
	private DFEVar getCurrentY(DFEVar row) {
		DFEType type = dfeInt(row.getType().getTotalBits() + 3);

		return 2 * (row.cast(type) - CoarseMeScalars.windowY/2);
	}

	/**
	 * This method finds the minimum SAD among the SADs calculated by the different pipes in the current cycle
	 * @param sads     a list of results of the outputs of computeSAD for all the pipes
	 * @param vectors  the motion vectors for the corresponding SADs
	 * @return         returns the minimum as a 3d array
	 */
	private DFEVar[][][] getCycleMinimum(List<DFEVar[][]> sads, DFEVar[][] vectors) {
		//first dim is sub-blocksize 4, 8, 16, 32, second dim is number of sub blocks of that size
		DFEVar[][][] output = new DFEVar[sads[0].length][sads[0][0].length][3];
		optimization.pushEnableBitGrowth(true);
		for (int i = 0; i < sads[0].length; i++) {
			int numSubBlocks = (m_scalars.blockSize * m_scalars.blockSize) / (1 << (2 * i + 4));
			for (int j = 0; j < numSubBlocks; j++) {
				output[i][j][0] = sads[0][i][j];
				output[i][j][1] = vectors[0][0];
				output[i][j][2] = vectors[0][1];
				for (int k = 1; k < sads.size(); k++) {
					DFEVar smaller = sads[k][i][j] < output[i][j][0];
					output[i][j][0] = smaller ? sads[k][i][j] : output[i][j][0];
					output[i][j][1] = smaller ? vectors[k][0] : output[i][j][1];
					output[i][j][2] = smaller ? vectors[k][1] : output[i][j][2];
				}
			}
		}
		optimization.popEnableBitGrowth();

		return output;
	}

	/**
	 * This methods compares the minimum in this cycle to the aggregate minimum of the search window and update the aggregate minimum.
	 * If this is the first cycle of processing a search window, the previous aggregate minimum is ignored
	 * @param startOfBlock  whether this is the first cycle of processing a new search window or not
	 * @param curr          current aggregate minimum
	 * @return              updated aggregate minimum
	 */
	private DFEVar[][][] getAggregateMinimum(DFEVar startOfBlock, DFEVar[][][] curr) {
		DFEVar[][][] output = new DFEVar[curr.length][curr[0].length][curr[0][0].length];

		for (int i = 0; i < curr.length; i++) {
			int numSubBlocks = (m_scalars.blockSize * m_scalars.blockSize) / (1 << (2 * i + 4));
			for (int j = 0; j < numSubBlocks; j++) {
				output[i][j][0] = curr[i][j][0].getType().newInstance(this);
				output[i][j][1] = curr[i][j][1].getType().newInstance(this);
				output[i][j][2] = curr[i][j][2].getType().newInstance(this);
				DFEVar prevCost = stream.offset(output[i][j][0], -1);
				DFEVar prevX    = stream.offset(output[i][j][1], -1);
				DFEVar prevY    = stream.offset(output[i][j][2], -1);

				optimization.pushPipeliningFactor(0.0);
				DFEVar takeNew = (curr[i][j][0] < prevCost) | startOfBlock;
				optimization.popPipeliningFactor();

				output[i][j][0] <== takeNew ? curr[i][j][0] : prevCost;
				output[i][j][1] <== takeNew ? curr[i][j][1] : prevX;
				output[i][j][2] <== takeNew ? curr[i][j][2] : prevY;
			}
		}

		return output;
	}

	/**
	 * Converts the minimum motion vectors to the output stream format
	 * @param minimum  minimum 3D array as produced by getMinimum() methods
	 * @return         an instance of MvType structure suitable to output the motion vectors to the stream
	 */
	@SuppressWarnings("unchecked")
	private DFEStruct getMotionVectors(DFEVar[][][] minimum) {
		DFEStruct output = getMvType(m_scalars.blockSize).newInstance(this);
		int index = 0;
		for (int i = 4; i <= m_scalars.blockSize; i *= 2) {
			int numSubBlocks = ((DFEVector<DFEVar>)output["best"+2*i]).getSize() / 2;
			for (int j = 0; j < numSubBlocks; j++) {
				((DFEVector<DFEVar>)output["best"+2*i])[2 * j]     <== minimum[index][j][1].cast(dfeInt(16));
				((DFEVector<DFEVar>)output["best"+2*i])[2 * j + 1] <== minimum[index][j][2].cast(dfeInt(16));
			}
			index++;
		}

		if (output.getType().getFieldNames().contains("padding")) {
			output["padding"] <== constant.var((DFEType) output["padding"].getType(), Bits.allZeros(output["padding"].getType().getTotalBits()));
		}

		return output;
	}

	/**
	 * Creates a custom structure type for the output motion vector stream
	 * @param blockSize   the length of a block side in pixels
	 * @return            DFE struct type that can contain all the motion vectors for the block and the sub-blocks
	 */
	public static DFEStructType getMvType(int blockSize) {
		List<StructFieldType> fields = new ArrayList<StructFieldType>();
		int totalBits = 0;
		for (int i = 4; i <= blockSize; i *= 2) {
			int numBlocks = (blockSize * blockSize) / (i * i);
			fields.add(sft("best" + i * 2, new DFEVectorType<DFEVar>(dfeInt(16), 2 * numBlocks)));//alternating X, Y
			totalBits += numBlocks * 32;
		}
		int padding = MathUtils.nextMultiple(totalBits, 128) - totalBits;
		if (padding > 0) {
			fields.add(sft("padding", dfeRawBits(padding)));
		}

		DFEStructType type = new DFEStructType(fields.toArray(new StructFieldType[1]));

		return type;
	}

}