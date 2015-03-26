=============================
MotionEstimation
=============================

This project implements a Full Search Block Matching Motion Estimation that can be used for video encoding.


Introduction
------------

Motion Estimation is used in video encoding to describe a video frame by motion vectors from other frames of the video. (Usually, adjacent frames) The motion vectors are defined for blocks of pixels. In this implementation, the algorithm can process 4x4, 8x8, and 16x16 blocks. The source frame and a single reference frames should be tiled as blocks before getting streamed into the DFE. The DFE performs a full search block matching motion estimation on all the blocks in a search window around the source block and outputs the best corresponding motion vector for that block.

Usage
-----

Import the MaxIDE project in the APP folder and compile and run either in simulation or on a DFE.

Features
--------

- This is a simple Full Search Block Matching Motion Estimation. No Intra-frame prediction technique is implemented.
- The algorithm expects halo blocks around the reference frame that ensure that the search windows will never go out of bound.
- The source and reference frames are two separate streams that must be provided by the CPU.
- The motion vectors are provided to the CPU as a stream.

