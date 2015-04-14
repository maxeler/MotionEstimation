package motionEstimation;


import java.util.ArrayList;
import java.util.List;

import com.maxeler.maxcompiler.v2.kernelcompiler.KernelLib;
import com.maxeler.maxcompiler.v2.kernelcompiler.RoundingMode;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.CounterChain;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.Mem.RamPortMode;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.Mem.RamPortParams;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.Mem.RamWriteMode;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVector;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVectorType;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

/**
 * This class abstracts the buffer that is used to store the active search window and provide the reference blocks to the Kernel.
 * Two public methods exists to interact with the buffer.
 *
 */
public class CoarseMeBuffer extends KernelLib {

	private final CoarseMeScalars m_scalars;
	private final DFEVar[][] m_bufferOutput;
	private final DFEVar m_outputEnable;
	private final DFEType m_addressType;

	/**
	 * This constructs the buffer. The buffer consists of FMem to hold the active search window, the address generators to read/write
	 * to the buffer, the controllers to disable/enable the buffer.
	 * @param scalars  The scalars variables object
	 */
	protected CoarseMeBuffer(CoarseMeScalars scalars) {
		super(scalars.getKernel());
		stream.suppressOffsetVectorWarnings();

		m_scalars = scalars;
		m_addressType = dfeUInt(MathUtils.bitsToAddress(scalars.bufferDepth));

		CoarseMeBufferControl buffCtrl = new CoarseMeBufferControl(scalars);
		m_outputEnable   = buffCtrl.getOutputEnable();

		DFEVectorType<DFEVar> type = new DFEVectorType<DFEVar>(dfeUInt(16), scalars.numPipes);
		DFEVector<DFEVar> reference = io.input("reference", type, buffCtrl.getReferenceInputEnable());

		DFEVar[] writeAddress = getWriteAddress(buffCtrl.getReferenceInputEnable());

		DFEVar[] readAddress1 = getReadAddresses(buffCtrl.getReadEnableBuffer1(), 0);
		DFEVar[] readAddress2 = getReadAddresses(buffCtrl.getReadEnableBuffer2(), 1);

		int cyclesToReadAhead = MathUtils.ceilDivide(scalars.blockSize + scalars.numPipes - 1, scalars.numPipes);

		DFEVar[][] bufferOutput = createBuffer(reference, buffCtrl.getReferenceInputEnable(),
				                               writeAddress, readAddress1, readAddress2,
				                               buffCtrl.getWhichBuffer(), cyclesToReadAhead);

		DFEVar rotate = buffCtrl.getWhichBuffer() ? stream.offset(readAddress2[scalars.blockSize], -cyclesToReadAhead)
                                               : stream.offset(readAddress1[scalars.blockSize], -cyclesToReadAhead);

		m_bufferOutput = rotateBlock(bufferOutput, rotate);


	}


	/**
	 * This method acts as the address generator for writing the reference input data into the buffer
	 * @param enable  The write enable signal for the buffer
	 * @return Returns the write address, and the row in which the data should be written to
	 */
	private DFEVar[] getWriteAddress(DFEVar enable) {
		CounterChain writeChain = control.count.makeCounterChain(enable);
		DFEVar wBlockRow = writeChain.addCounter(m_scalars.heightPlusHaloInBlocks, 1);
		DFEVar wBlockCol = writeChain.addCounter(m_scalars.widthPlusHaloInBlocks, 1);
		DFEVar wRow      = writeChain.addCounter(m_scalars.blockSize, 1);
		DFEVar wColumn   = m_scalars.blockSize > m_scalars.numPipes ? writeChain.addCounter(m_scalars.blockSize / m_scalars.numPipes, 1) : constant.var(0);

		/* addresses are in sequence with the order that data arrived.
		 * The address is going to wrap back to 0 when it reaches the maximum address in the buffer.
		 * By that time, the previous data in the buffer is no longer useful */
		optimization.pushEnableBitGrowth(true);
		optimization.pushRoundingMode(RoundingMode.TRUNCATE);
		DFEVar writeAddress = (wColumn + wBlockCol * (m_scalars.blockSize / m_scalars.numPipes)
		                               + wBlockRow * m_scalars.widthPlusHaloInBlocks * (m_scalars.blockSize / m_scalars.numPipes)).cast(m_addressType);
		optimization.popRoundingMode();
		optimization.popEnableBitGrowth();

		DFEVar[] output = new DFEVar[2];
		output[0] = writeAddress;
		output[1] = wRow;

		return output;
	}


	/**
	 * Creates the buffer using FMem. The memory organised as an array of memories for each row in a block to allow the kernel
	 * to read a column of a block in one cycle. There are two reading ports to allow the blocks to be read from two different
	 * rows at the same time.
	 * @param reference           The current reference input received from CPU
	 * @param writeEnable         Write enable signal
	 * @param writeAddress        Write address generated by getWriteAddress() method
	 * @param readAddress1        Read address generated by getReadAddresses() method
	 * @param readAddress2        Read address generated by getReadAddresses() method
	 * @param whichBuffer         Selects which of the two buffer to read from
	 * @param cyclesToReadAhead   The number of cycles that it takes to read the columns from memory to form the blocks for all the pipes
	 * @return                    The output is a 2d array, wider than block size, to contain all the blocks for different pipes
	 */
	private DFEVar[][] createBuffer(DFEVector<DFEVar> reference,
			                        DFEVar writeEnable,
			                        DFEVar[] writeAddress,
			                        DFEVar[] readAddress1,
			                        DFEVar[] readAddress2,
			                        DFEVar whichBuffer,
			                        int cyclesToReadAhead)
	{
		DFEVar[][] bufferOutput = new DFEVar[m_scalars.blockSize][m_scalars.blockSize + m_scalars.numPipes - 1];
		for (int i = 0; i < m_scalars.blockSize; i++) {
			RamPortParams<DFEVector<DFEVar>> writePort = mem.makeRamPortParams(RamPortMode.WRITE_ONLY, writeAddress[0], reference.getType())
			                                                .withDataIn(reference)
                                                            .withWriteEnable((writeAddress[1] === i) & writeEnable);

			RamPortParams<DFEVector<DFEVar>> readPort1  = mem.makeRamPortParams(RamPortMode.READ_ONLY, readAddress1[i], reference.getType());
			RamPortParams<DFEVector<DFEVar>> readPort2  = mem.makeRamPortParams(RamPortMode.READ_ONLY, readAddress2[i], reference.getType());

			DFEVector<DFEVar> ram1 = mem.ramDualPort(m_scalars.bufferDepth, RamWriteMode.READ_FIRST, writePort, readPort1).getOutputB();
			DFEVector<DFEVar> ram2 = mem.ramDualPort(m_scalars.bufferDepth, RamWriteMode.READ_FIRST, writePort, readPort2).getOutputB();

			for (int j = 0; j < cyclesToReadAhead; j++) {
				DFEVector<DFEVar> temp = whichBuffer ? stream.offset(ram2, j - cyclesToReadAhead + 1) : stream.offset(ram1, j - cyclesToReadAhead + 1);
				for (int k = 0; k < m_scalars.numPipes; k++) {
					if (j * m_scalars.numPipes + k < bufferOutput[0].length) {
						bufferOutput[i][j * m_scalars.numPipes + k] = temp[k];
					}
				}
			}
		}

		return bufferOutput;
	}

	/**
	 * This method rotates the rows in the output block of the buffer, this is needed because of the way the memory is organised
	 * @param input   output of the buffer
	 * @param rotate  amount of rotation
	 * @return        rotated block
	 */
	private DFEVar[][] rotateBlock(DFEVar[][] input, DFEVar rotate) {
		DFEVectorType<DFEVar> colType = new DFEVectorType<DFEVar>(input[0][0].getType(), input.length);
		DFEVar[][] output = new DFEVar[input.length][input[0].length];
		for (int i = 0; i < input[0].length; i++) {
			List<DFEVar> temp = new ArrayList<DFEVar>(input.length);
			for (int j = 0; j < input.length; j++) {
				temp.add(input[j][i]);
			}

			DFEVector<DFEVar> column = colType.newInstance(this, temp);
			column = column.rotateElementsRight(rotate);
			for (int j = 0; j < input.length; j++) {
				output[j][i] = column[j];
			}
		}

		return output;
	}

	/**
	 * This method acts as the address generator for writing the reference input data into the buffer
	 * @param enable  read enable signal
	 * @param buffer  determines which buffer to read from, can be 0 or 1
	 * @return        returns the read address for the current cycle
	 */
	private DFEVar[] getReadAddresses(DFEVar enable, int buffer) {
		int totalWindowWidth = MathUtils.ceilDivide(CoarseMeScalars.windowX + m_scalars.blockSize - 1, m_scalars.numPipes);

		CounterChain readChain = control.count.makeCounterChain(enable);
		DFEVar rBlockRow = readChain.addCounter(m_scalars.heightInBlocks, 1);
		DFEVar rBlockCol = readChain.addCounter(m_scalars.widthInBlocks, 1);
		DFEVar rRow      = (readChain.addCounter(CoarseMeScalars.windowY, 2) + buffer);
		DFEVar rColumn   = readChain.addCounter(totalWindowWidth, 1);

		optimization.pushEnableBitGrowth(true);
		optimization.pushRoundingMode(RoundingMode.TRUNCATE);
		rRow = rRow + (m_scalars.haloBlocksY * m_scalars.blockSize) - (CoarseMeScalars.windowY/2);
		rColumn = rColumn + (m_scalars.haloBlocksX * m_scalars.blockSize) - (CoarseMeScalars.windowX/2);
		DFEVar rotate = rRow.cast(dfeUInt(MathUtils.bitsToAddress(m_scalars.blockSize)));

		DFEVar rY = ((rBlockRow * m_scalars.blockSize + rRow) / m_scalars.blockSize).cast(m_scalars.heightInBlocks.getType());
		DFEVar rX = rColumn + rBlockCol * (m_scalars.blockSize / m_scalars.numPipes);

		DFEVar readAddress1 = (rX + rY       * m_scalars.widthPlusHaloInBlocks * (m_scalars.blockSize / m_scalars.numPipes)).cast(m_addressType);
		DFEVar readAddress2 = (rX + (rY + 1) * m_scalars.widthPlusHaloInBlocks * (m_scalars.blockSize / m_scalars.numPipes)).cast(m_addressType);
		optimization.popRoundingMode();
		optimization.popEnableBitGrowth();

		DFEVar[] readAddresses = new DFEVar[m_scalars.blockSize + 1];
		for (int i = 0; i < m_scalars.blockSize; i++) {
			readAddresses[i] = i < rotate ? readAddress2 : readAddress1;
		}

		readAddresses[m_scalars.blockSize] = rotate;


		return readAddresses;
	}


	/**
	 * @param pipe  the pipe that processes the block
	 * @return      current reference block for the corresponding pipe from the output of the buffer
	 */
	public DFEVar[][] getReferenceBlock(int pipe) {
		DFEVar[][] output = new DFEVar[m_scalars.blockSize][m_scalars.blockSize];
		for (int i = 0; i < m_scalars.blockSize; i++) {
			for (int j = 0; j < m_scalars.blockSize; j++) {
				output[i][j] = m_bufferOutput[i][j + pipe];
			}
		}
		return output;
	}

	/**
	 * @return  the output enable signal
	 */
	public DFEVar getOutputEnable() {
		return m_outputEnable;
	}

}