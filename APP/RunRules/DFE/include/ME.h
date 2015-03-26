/**\file */
#ifndef SLIC_DECLARATIONS_ME_H
#define SLIC_DECLARATIONS_ME_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define ME_blockSize (16)
#define ME_numPipes (2)
#define ME_windowX (64)
#define ME_windowY (32)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] ticks_CoarseMeKernel The number of ticks for which kernel "CoarseMeKernel" will run.
 * \param [in] inscalar_CoarseMeKernel_blocksToComplete Input scalar parameter "CoarseMeKernel.blocksToComplete".
 * \param [in] inscalar_CoarseMeKernel_fillBlocks Input scalar parameter "CoarseMeKernel.fillBlocks".
 * \param [in] inscalar_CoarseMeKernel_heightInBlocks Input scalar parameter "CoarseMeKernel.heightInBlocks".
 * \param [in] inscalar_CoarseMeKernel_lambda Input scalar parameter "CoarseMeKernel.lambda".
 * \param [in] inscalar_CoarseMeKernel_numBlocks Input scalar parameter "CoarseMeKernel.numBlocks".
 * \param [in] inscalar_CoarseMeKernel_pictureHeight Input scalar parameter "CoarseMeKernel.pictureHeight".
 * \param [in] inscalar_CoarseMeKernel_widthInBlocks Input scalar parameter "CoarseMeKernel.widthInBlocks".
 * \param [in] inscalar_CoarseMeKernel_widthPlusHaloInBlocks Input scalar parameter "CoarseMeKernel.widthPlusHaloInBlocks".
 * \param [in] instream_reference Stream "reference".
 * \param [in] instream_size_reference The size of the stream instream_reference in bytes.
 * \param [in] instream_source Stream "source".
 * \param [in] instream_size_source The size of the stream instream_source in bytes.
 * \param [out] outstream_mv Stream "mv".
 * \param [in] outstream_size_mv The size of the stream outstream_mv in bytes.
 */
void ME(
	uint64_t ticks_CoarseMeKernel,
	uint64_t inscalar_CoarseMeKernel_blocksToComplete,
	uint64_t inscalar_CoarseMeKernel_fillBlocks,
	uint64_t inscalar_CoarseMeKernel_heightInBlocks,
	uint64_t inscalar_CoarseMeKernel_lambda,
	uint64_t inscalar_CoarseMeKernel_numBlocks,
	uint64_t inscalar_CoarseMeKernel_pictureHeight,
	uint64_t inscalar_CoarseMeKernel_widthInBlocks,
	uint64_t inscalar_CoarseMeKernel_widthPlusHaloInBlocks,
	const void *instream_reference,
	size_t instream_size_reference,
	const void *instream_source,
	size_t instream_size_source,
	void *outstream_mv,
	size_t outstream_size_mv);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] ticks_CoarseMeKernel The number of ticks for which kernel "CoarseMeKernel" will run.
 * \param [in] inscalar_CoarseMeKernel_blocksToComplete Input scalar parameter "CoarseMeKernel.blocksToComplete".
 * \param [in] inscalar_CoarseMeKernel_fillBlocks Input scalar parameter "CoarseMeKernel.fillBlocks".
 * \param [in] inscalar_CoarseMeKernel_heightInBlocks Input scalar parameter "CoarseMeKernel.heightInBlocks".
 * \param [in] inscalar_CoarseMeKernel_lambda Input scalar parameter "CoarseMeKernel.lambda".
 * \param [in] inscalar_CoarseMeKernel_numBlocks Input scalar parameter "CoarseMeKernel.numBlocks".
 * \param [in] inscalar_CoarseMeKernel_pictureHeight Input scalar parameter "CoarseMeKernel.pictureHeight".
 * \param [in] inscalar_CoarseMeKernel_widthInBlocks Input scalar parameter "CoarseMeKernel.widthInBlocks".
 * \param [in] inscalar_CoarseMeKernel_widthPlusHaloInBlocks Input scalar parameter "CoarseMeKernel.widthPlusHaloInBlocks".
 * \param [in] instream_reference Stream "reference".
 * \param [in] instream_size_reference The size of the stream instream_reference in bytes.
 * \param [in] instream_source Stream "source".
 * \param [in] instream_size_source The size of the stream instream_source in bytes.
 * \param [out] outstream_mv Stream "mv".
 * \param [in] outstream_size_mv The size of the stream outstream_mv in bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *ME_nonblock(
	uint64_t ticks_CoarseMeKernel,
	uint64_t inscalar_CoarseMeKernel_blocksToComplete,
	uint64_t inscalar_CoarseMeKernel_fillBlocks,
	uint64_t inscalar_CoarseMeKernel_heightInBlocks,
	uint64_t inscalar_CoarseMeKernel_lambda,
	uint64_t inscalar_CoarseMeKernel_numBlocks,
	uint64_t inscalar_CoarseMeKernel_pictureHeight,
	uint64_t inscalar_CoarseMeKernel_widthInBlocks,
	uint64_t inscalar_CoarseMeKernel_widthPlusHaloInBlocks,
	const void *instream_reference,
	size_t instream_size_reference,
	const void *instream_source,
	size_t instream_size_source,
	void *outstream_mv,
	size_t outstream_size_mv);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t ticks_CoarseMeKernel; /**<  [in] The number of ticks for which kernel "CoarseMeKernel" will run. */
	uint64_t inscalar_CoarseMeKernel_blocksToComplete; /**<  [in] Input scalar parameter "CoarseMeKernel.blocksToComplete". */
	uint64_t inscalar_CoarseMeKernel_fillBlocks; /**<  [in] Input scalar parameter "CoarseMeKernel.fillBlocks". */
	uint64_t inscalar_CoarseMeKernel_heightInBlocks; /**<  [in] Input scalar parameter "CoarseMeKernel.heightInBlocks". */
	uint64_t inscalar_CoarseMeKernel_lambda; /**<  [in] Input scalar parameter "CoarseMeKernel.lambda". */
	uint64_t inscalar_CoarseMeKernel_numBlocks; /**<  [in] Input scalar parameter "CoarseMeKernel.numBlocks". */
	uint64_t inscalar_CoarseMeKernel_pictureHeight; /**<  [in] Input scalar parameter "CoarseMeKernel.pictureHeight". */
	uint64_t inscalar_CoarseMeKernel_widthInBlocks; /**<  [in] Input scalar parameter "CoarseMeKernel.widthInBlocks". */
	uint64_t inscalar_CoarseMeKernel_widthPlusHaloInBlocks; /**<  [in] Input scalar parameter "CoarseMeKernel.widthPlusHaloInBlocks". */
	const void *instream_reference; /**<  [in] Stream "reference". */
	size_t instream_size_reference; /**<  [in] The size of the stream instream_reference in bytes. */
	const void *instream_source; /**<  [in] Stream "source". */
	size_t instream_size_source; /**<  [in] The size of the stream instream_source in bytes. */
	void *outstream_mv; /**<  [out] Stream "mv". */
	size_t outstream_size_mv; /**<  [in] The size of the stream outstream_mv in bytes. */
} ME_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void ME_run(
	max_engine_t *engine,
	ME_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *ME_run_nonblock(
	max_engine_t *engine,
	ME_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void ME_run_group(max_group_t *group, ME_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *ME_run_group_nonblock(max_group_t *group, ME_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void ME_run_array(max_engarray_t *engarray, ME_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *ME_run_array_nonblock(max_engarray_t *engarray, ME_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* ME_convert(max_file_t *maxfile, ME_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* ME_init(void);

/* Error handling functions */
int ME_has_errors(void);
const char* ME_get_errors(void);
void ME_clear_errors(void);
/* Free statically allocated maxfile data */
void ME_free(void);
/* These are dummy functions for hardware builds. */
int ME_simulator_start(void);
int ME_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_ME_H */

