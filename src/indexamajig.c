/*
 * indexamajigFork.c
 *
 * Index patterns, output hkl+intensity etc.
 *
 * (c) 2006-2012 Thomas White <taw@physics.org>
 * (c) 2012-     Chunhong Yoon <chun.hong.yoon@desy.de>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <hdf5.h>
#include <gsl/gsl_errno.h>
#include <pthread.h>
#include <sys/wait.h>
#include <fcntl.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <crystfel/utils.h>
#include <crystfel/hdf5-file.h>
#include <crystfel/index.h>
#include <crystfel/peaks.h>
#include <crystfel/detector.h>
#include <crystfel/filters.h>
#include <crystfel/thread-pool.h>
#include <crystfel/beam-parameters.h>
#include <crystfel/geometry.h>
#include <crystfel/stream.h>
#include <crystfel/reflist-utils.h>


/* Write statistics at APPROXIMATELY this interval */
#define STATS_EVERY_N_SECONDS (2)

#define LINE_LENGTH 1024

#define BUFSIZE 64

#ifdef HAVE_CLOCK_GETTIME
static double get_time()
{
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);
	double sec = (double) tp.tv_sec+ (double) tp.tv_nsec/1000000000;
	return sec; //nano resolution
}
#else
/* Fallback version of the above.  The time according to gettimeofday() is not
 * monotonic, so measuring intervals based on it will screw up if there's a
 * timezone change (e.g. daylight savings) while the program is running. */
static double get_time()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	double sec = (double) tp.tv_sec+ (double) tp.tv_usec/1000000;
	return sec; //micro resolution
}
#endif

enum {
	PEAK_ZAEF,
	PEAK_HDF5,
};


/* Information about the indexing process which is common to all patterns */
struct static_index_args
{
	UnitCell *cell;
	int config_cmfilter;
	int config_noisefilter;
	int config_verbose;
	int stream_flags;         /* What goes into the output? */
	int config_polar;
	int config_satcorr;
	int config_closer;
	int config_insane;
	int config_bgsub;
	float threshold;
	float min_gradient;
	float min_snr;
	double min_int_snr;
	struct detector *det;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	int peaks;                /* Peak detection method */
	int cellr;
	float tols[4];
	struct beam_params *beam;
	const char *element;
	const char *hdf5_peak_path;
	double ir_inn;
	double ir_mid;
	double ir_out;
	
	/* Output stream */
	pthread_mutex_t *output_mutex;  /* Protects the output stream */
	FILE *ofh;

	char *outfile;	
	const struct copy_hdf5_field *copyme;
	int nProcesses;
};


/* Information about the indexing process for one pattern */
struct index_args
{
	/* "Input" */
	char *filename;
	struct static_index_args static_args;

	/* "Output" */
	int indexable;
};


/* Information needed to choose the next task and dispatch it */
struct queue_args
{
	FILE *fh;
	char *prefix;
	int config_basename;
	struct static_index_args static_args;

	char *use_this_one_instead;

	int n_indexable;
	int n_processed;
	int n_indexable_last_stats;
	int n_processed_last_stats;

	int n_indexableTotal;
	int n_processedTotal;
	int n_indexable_last_statsTotal;
	int n_processed_last_statsTotal;

	int nPerProcess;
	int updateReader;

	double t_last_stats;

};

// Count number of patterns in .lst
int count_patterns(FILE *fh){
    char *rval;
    int n_total = 0;
    do{
        char line[LINE_LENGTH];
        rval = fgets(line,LINE_LENGTH-1,fh);
        if (rval != NULL){
            n_total++;
        }
    }while(rval!=NULL);
    if (ferror(fh)) {
        printf("Read error\n");
        return -1;
    }
    return n_total;
}

// Assign a batch number to a process
int getBatchNum(int pid[], int length){
    int i, id = 0;
    for (i=0; i<length; i++) {
        id += ((pid[i]>0)*pow(2,i)); // child = 1, parent = 0
    }
    return id;
}

static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Process and index FEL diffraction images.\n"
"\n"
" -h, --help               Display this help message.\n"
"\n"
" -i, --input=<filename>   Specify file containing list of images to process.\n"
"                           '-' means stdin, which is the default.\n"
" -o, --output=<filename>  Write output stream to this file. '-' for stdout.\n"
"                           Default: indexamajig.stream\n"
"\n"
"     --indexing=<methods> Use 'methods' for indexing.  Provide one or more\n"
"                           methods separated by commas.  Choose from:\n"
"                            none     : no indexing (default)\n"
"                            dirax    : invoke DirAx\n"
"                            mosflm   : invoke MOSFLM (DPS)\n"
"                            reax     : DPS algorithm with known unit cell\n"
" -g. --geometry=<file>    Get detector geometry from file.\n"
" -b, --beam=<file>        Get beam parameters from file (provides nominal\n"
"                           wavelength value if no per-shot value is found in\n"
"                           the HDF5 files.\n"
" -p, --pdb=<file>         PDB file from which to get the unit cell to match.\n"
"                           Default: 'molecule.pdb'.\n"
"     --basename           Remove the directory parts of the filenames.\n"
" -x, --prefix=<p>         Prefix filenames from input file with <p>.\n"
"     --peaks=<method>     Use 'method' for finding peaks.  Choose from:\n"
"                           zaef  : Use Zaefferer (2000) gradient detection.\n"
"                                    This is the default method.\n"
"                           hdf5  : Get from a table in HDF5 file.\n"
"     --hdf5-peaks=<p>     Find peaks table in HDF5 file here.\n"
"                           Default: /processing/hitfinder/peakinfo\n"
"\n\n"
"You can control what information is included in the output stream using\n"
"' --record=<flag1>,<flag2>,<flag3>' and so on.  Possible flags are:\n\n"
" integrated        Include a list of reflection intensities, produced by\n"
"                    integrating around predicted peak locations.\n"
"\n"
" peaks             Include peak locations and intensities from the peak\n"
"                    search.\n"
"\n"
" peaksifindexed    As 'peaks', but only if the pattern could be indexed.\n"
"\n"
" peaksifnotindexed As 'peaks', but only if the pattern could NOT be indexed.\n"
"\n\n"
"The default is '--record=integrated'.\n"
"\n\n"
"For more control over the process, you might need:\n\n"
"     --cell-reduction=<m> Use <m> as the cell reduction method. Choose from:\n"
"                           none    : no matching, just use the raw cell.\n"
"                           reduce  : full cell reduction.\n"
"                           compare : match by at most changing the order of\n"
"                                     the indices.\n"
"                           compare_ab : compare 'a' and 'b' lengths only.\n"
"     --tolerance=<a,b,c,angl>  Set the tolerance for a,b,c axis (in %%)\n"
"                                and for the angles (in deg) when reducing\n"
"                                or comparing (default is 5%% and 1.5deg)\n"
"     --filter-cm          Perform common-mode noise subtraction on images\n"
"                           before proceeding.  Intensities will be extracted\n"
"                           from the image as it is after this processing.\n"
"     --filter-noise       Apply an aggressive noise filter which sets all\n"
"                           pixels in each 3x3 region to zero if any of them\n"
"                           have negative values.  Intensity measurement will\n"
"                           be performed on the image as it was before this.\n"
"     --unpolarized        Don't correct for the polarisation of the X-rays.\n"
"     --no-sat-corr        Don't correct values of saturated peaks using a\n"
"                           table included in the HDF5 file.\n"
"     --threshold=<n>      Only accept peaks above <n> ADU.  Default: 800.\n"
"     --min-gradient=<n>   Minimum gradient for Zaefferer peak search.\n"
"                           Default: 100,000.\n"
"     --min-snr=<n>        Minimum signal-to-noise ratio for peaks.\n"
"                           Default: 5.\n"
"     --min-integration-snr=<n>  Minimum signal-to-noise ratio for peaks\n"
"                                 during integration. Default: -infinity.\n"
" -e, --image=<element>    Use this image from the HDF5 file.\n"
"                           Example: /data/data0.\n"
"                           Default: The first one found.\n"
"\n"
"\nFor time-resolved stuff, you might want to use:\n\n"
"     --copy-hdf5-field <f>  Copy the value of field <f> into the stream. You\n"
"                             can use this option as many times as you need.\n"
"\n"
"\nOptions for greater performance or verbosity:\n\n"
"     --verbose            Be verbose about indexing.\n"
" -j <n>                   Run <n> analyses in parallel.  Default 1.\n"
"\n"
"\nOptions you probably won't need:\n\n"
"     --no-check-prefix    Don't attempt to correct the --prefix.\n"
"     --no-closer-peak     Don't integrate from the location of a nearby peak\n"
"                           instead of the position closest to the reciprocal\n"
"                           lattice point.\n"
"     --insane             Don't check that the reduced cell accounts for at\n"
"                           least 10%% of the located peaks.\n"
"     --no-bg-sub          Don't subtract local background estimates from\n"
"                           integrated intensities.\n"
"\n"
"\nYou can tune the CPU affinities for enhanced performance on NUMA machines:\n"
"\n"
"     --cpus=<n>           Specify number of CPUs.  This is NOT the same as\n"
"                           giving the number of analyses to run in parallel.\n"
"     --cpugroup=<n>       Batch threads in groups of this size.\n"
"     --cpuoffset=<n>      Start using CPUs at this group number.\n"
);
}

int readUpdate(void *qp, int fd_pipe[][2], int finish){
	struct queue_args *qargs = qp;
	int i,j;
	int rFlag;
	int pipesOpen = 0;
	char bufferR[BUFSIZE];			
	memset(bufferR, 0, BUFSIZE);
	int numFields = 2;

	for (i=0; i<qargs->static_args.nProcesses; i++){
		if (finish && i == qargs->updateReader) {
			// do nothing
		}else{
			rFlag = read(fd_pipe[i][0], bufferR, BUFSIZE-1); // wait till something to read
			if (rFlag != 0){
				char * pch; 
				char delims[] = "#";
  				pch = strtok (bufferR, delims); 
  				for (j=0;j<numFields;j++) 
  				{ 
					switch (j){
						case 0:
							qargs->n_indexableTotal += atoi(pch);	
							break;
						case 1:
							qargs->n_processedTotal += atoi(pch);	
						break;
					}
    				pch = strtok (NULL, delims); 
  				} 
				memset(bufferR, 0, BUFSIZE);
				pipesOpen++;
			}
		}
	}
	STATUS("%i out of %i indexed so far,"
		       " %i out of %i since the last message.\n\n",
			qargs->n_indexableTotal, qargs->n_processedTotal,
			qargs->n_indexableTotal-qargs->n_indexable_last_statsTotal,
			qargs->n_processedTotal-qargs->n_processed_last_statsTotal);

	qargs->n_indexable_last_statsTotal = qargs->n_indexableTotal;
	qargs->n_processed_last_statsTotal = qargs->n_processedTotal;

	return pipesOpen;	
}

// Manipulate image
static void process_image(char ***array, int batch, void *qp, int fd_pipe[][2]) {
	
	struct queue_args *qargs = qp;
	float *data_for_measurement;
	size_t data_size;
	UnitCell *cell = qargs->static_args.cell;
	int config_cmfilter = qargs->static_args.config_cmfilter;
	int config_noisefilter = qargs->static_args.config_noisefilter;
	int config_verbose = qargs->static_args.config_verbose;
	//int config_polar = qargs->static_args.config_polar;
	IndexingMethod *indm = qargs->static_args.indm;
	struct beam_params *beam = qargs->static_args.beam;

	int i;	
	int r, check;
	struct hdfile *hdfile;
	char *outfile = qargs->static_args.outfile;
	int nPerBatch = qargs->nPerProcess;

	for (i=0; i<nPerBatch; i++) {
		struct image image;
		image.features = NULL;
		image.data = NULL;
		image.flags = NULL;
		image.indexed_cell = NULL;
		image.det = copy_geom(qargs->static_args.det);
		image.copyme = qargs->static_args.copyme;
		image.beam = beam;
		image.id = batch;	// MUST SET ID FOR MOSFLM TO WORK PROPERLY

		if ( beam == NULL ) {
			ERROR("Warning: no beam parameters file.\n");
			ERROR("I'm going to assume 1 ADU per photon, which is almost");
			ERROR(" certainly wrong.  Peak sigmas will be incorrect.\n");
		}

		//pargs->indexable = 0;

		char *filename = NULL;
		char *line = array[batch][i];
		chomp(line);
//printf("%d ***************%s\n",batch,line);
		filename = malloc(strlen(qargs->prefix)+strlen(line)+1);
		snprintf(filename,LINE_LENGTH-1,"%s%s",qargs->prefix,line); //maximum print length
		free(line);
		image.filename = filename;
//printf("%d ***************%s\n",batch,filename);
		hdfile = hdfile_open(filename);
		if (hdfile == NULL) return;

		r = hdfile_set_first_image(hdfile, "/"); // Need this to read hdf5 files
		if (r) {
			ERROR("Couldn't select first path\n");
			hdfile_close(hdfile);
			return;
		}

		check = hdf5_read(hdfile,&image,1);
		if (check == 1){
			hdfile_close(hdfile);
			return;
		}

		if ( (image.width != image.det->max_fs+1)
		  || (image.height != image.det->max_ss+1) )
		{
			ERROR("Image size doesn't match geometry size"
			      " - rejecting image.\n");
			ERROR("Image size: %i,%i.  Geometry size: %i,%i\n",
			      image.width, image.height,
			      image.det->max_fs+1, image.det->max_ss+1);
			hdfile_close(hdfile);
			free_detector_geometry(image.det);
			return;
		}

		if ( image.lambda < 0.0 ) {
			if ( beam != NULL ) {
				ERROR("Using nominal photon energy of %.2f eV\n",
    	                          beam->photon_energy);
				image.lambda = ph_en_to_lambda(
			                          eV_to_J(beam->photon_energy));
			} else {
				ERROR("No wavelength in file, so you need to give "
			      "a beam parameters file with -b.\n");
				hdfile_close(hdfile);
				free_detector_geometry(image.det);
				return;
			}
		}
		fill_in_values(image.det, hdfile);

		if ( config_cmfilter ) {
			filter_cm(&image);
		}

		/* Take snapshot of image after CM subtraction but before
		 * the aggressive noise filter. */
		data_size = image.width*image.height*sizeof(float);
		data_for_measurement = malloc(data_size);

		if ( config_noisefilter ) {
			filter_noise(&image, data_for_measurement);
		} else {
			memcpy(data_for_measurement, image.data, data_size);
		}

		switch ( qargs->static_args.peaks )
		{
		case PEAK_HDF5 :
			/* Get peaks from HDF5 */
			if ( get_peaks(&image, hdfile,
		               qargs->static_args.hdf5_peak_path) )
			{
				ERROR("Failed to get peaks from HDF5 file.\n");
			}
			break;
		case PEAK_ZAEF :
			search_peaks(&image, qargs->static_args.threshold,
		             qargs->static_args.min_gradient,
		             qargs->static_args.min_snr,
		             qargs->static_args.ir_inn,
		             qargs->static_args.ir_mid,
		             qargs->static_args.ir_out);
			break;
		}
		
		/* Get rid of noise-filtered version at this point
		 * - it was strictly for the purposes of peak detection. */
		free(image.data);
		image.data = data_for_measurement;

		/* Calculate orientation matrix (by magic) */
		image.div = beam->divergence;
		image.bw = beam->bandwidth;
		image.profile_radius = 0.0001e9;

		///// RUN INDEXING HERE /////
		index_pattern(&image, cell, indm, qargs->static_args.cellr,
		     config_verbose, qargs->static_args.ipriv,
		     qargs->static_args.config_insane, qargs->static_args.tols);

		if ( image.indexed_cell != NULL ) {
			//pargs->indexable = 1;
			image.reflections = find_intersections(&image,
		                                       image.indexed_cell);
			if ( image.reflections != NULL ) {
				integrate_reflections(&image,
					      qargs->static_args.config_closer,
					      qargs->static_args.config_bgsub,
					      qargs->static_args.min_int_snr,
					      qargs->static_args.ir_inn,
					      qargs->static_args.ir_mid,
					      qargs->static_args.ir_out);
			}
		} else {
			image.reflections = NULL;
		}

		// Write Lock
		struct flock fl = {F_WRLCK, SEEK_SET,   0,      0,     0 };
    	int fd;
	    fl.l_pid = getpid();
	
		char *outfilename = NULL;
		chomp(outfile); // may not need this
//printf("%d ***************%s\n",batch,outfile);		
		outfilename = malloc(strlen(outfile)+1);
		snprintf(outfilename,LINE_LENGTH-1,"%s",outfile); //maximum print length
//printf("%d ***************%s\n",batch,outfilename);
   		if ((fd = open(outfilename, O_WRONLY)) == -1) {
        	perror("Error on opening\n");
        	exit(1);
    	}

	    if (fcntl(fd, F_SETLKW, &fl) == -1) {
	        perror("Error on setting lock wait\n");
	        exit(1);
	    }

		// LOCKED! Write chunk
		FILE *fh;
//printf("%d ***************%s\n",batch,outfilename);
		fh = fopen(outfilename,"a");
		if (fh == NULL) { 
			perror("Error inside lock\n"); 
		}
		write_chunk(fh, &image, hdfile, qargs->static_args.stream_flags);		
		fclose(fh);

		// Unlock stream for other processes
	    fl.l_type = F_UNLCK;  // set to unlock same region
  	  	if (fcntl(fd, F_SETLK, &fl) == -1) {
   	    	perror("fcntl");
   	     	exit(1);
    	}
    	close(fd);

		///// WRITE UPDATE /////
		double seconds;
		qargs->n_indexable += ( image.indexed_cell != NULL );
		qargs->n_processed++;
		seconds = get_time();
		if ( seconds >= qargs->t_last_stats+STATS_EVERY_N_SECONDS 
			|| qargs->n_processed == qargs->nPerProcess) { // Must write if finished

			// WRITE PIPE HERE
			char bufferW[BUFSIZE];
			memset(bufferW, 0, BUFSIZE);
			sprintf(bufferW,"%d#%d",
				qargs->n_indexable - qargs->n_indexable_last_stats, 
				qargs->n_processed - qargs->n_processed_last_stats);
			write(fd_pipe[batch][1], bufferW, BUFSIZE-1);

			// Update stats
			qargs->n_processed_last_stats = qargs->n_processed;
			qargs->n_indexable_last_stats = qargs->n_indexable;
			qargs->t_last_stats = seconds;

			///// READ UPDATE /////
			if (batch == qargs->updateReader) {
				readUpdate(qargs,fd_pipe, 0);
			}
		}

		///// FREE /////
		/* Only free cell if found */
		cell_free(image.indexed_cell);
		reflist_free(image.reflections);
		free(image.data);
		if ( image.flags != NULL ) free(image.flags);
		image_feature_list_free(image.features);
		hdfile_close(hdfile);
		free_detector_geometry(image.det);
	}
}

static int parse_cell_reduction(const char *scellr, int *err,
                                int *reduction_needs_cell)
{
	*err = 0;
	if ( strcmp(scellr, "none") == 0 ) {
		*reduction_needs_cell = 0;
		return CELLR_NONE;
	} else if ( strcmp(scellr, "reduce") == 0) {
		*reduction_needs_cell = 1;
		return CELLR_REDUCE;
	} else if ( strcmp(scellr, "compare") == 0) {
		*reduction_needs_cell = 1;
		return CELLR_COMPARE;
	} else if ( strcmp(scellr, "compare_ab") == 0) {
		*reduction_needs_cell = 1;
		return CELLR_COMPARE_AB;
	} else {
		*err = 1;
		*reduction_needs_cell = 0;
		return CELLR_NONE;
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *outfile = NULL;
	FILE *fh;
	FILE *ofh;
	char *rval = NULL;
	int n_images;
	int config_noindex = 0;
	int config_cmfilter = 0;
	int config_noisefilter = 0;
	int config_verbose = 0;
	int config_polar = 1;
	int config_satcorr = 1;
	int config_checkprefix = 1;
	int config_closer = 1;
	int config_insane = 0;
	int config_bgsub = 1;
	int config_basename = 0;
	float threshold = 800.0;
	float min_gradient = 100000.0;
	float min_snr = 5.0;
	double min_int_snr = -INFINITY;
	struct detector *det;
	char *geometry = NULL;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	int indexer_needs_cell;
	int reduction_needs_cell;
	char *indm_str = NULL;
	UnitCell *cell;
	char *pdb = NULL;
	char *prefix = NULL;
	char *speaks = NULL;
	char *scellr = NULL;
	char *toler = NULL;
	float tols[4] = {5.0, 5.0, 5.0, 1.5}; /* a,b,c,angles (%,%,%,deg) */
	int cellr;
	int peaks;
	int nProcesses = 1;
	pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
	char *prepare_line;
	char prepare_filename[1024];
	struct queue_args qargs;
	struct beam_params *beam = NULL;
	char *element = NULL;
	double nominal_photon_energy;
	int stream_flags = STREAM_INTEGRATED;
	int cpu_num = 0;
	int cpu_groupsize = 1;
	int cpu_offset = 0;
	char *endptr;
	char *hdf5_peak_path = NULL;
	struct copy_hdf5_field *copyme;
	char *intrad = NULL;
	float ir_inn = 4.0;
	float ir_mid = 5.0;
	float ir_out = 7.0;
	
	copyme = new_copy_hdf5_field_list();
	if ( copyme == NULL ) {
		ERROR("Couldn't allocate HDF5 field list.\n");
		return 1;
	}

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"no-index",           0, &config_noindex,     1},
		{"indexing",           1, NULL,               'z'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"filter-cm",          0, &config_cmfilter,    1},
		{"filter-noise",       0, &config_noisefilter, 1},
		{"verbose",            0, &config_verbose,     1},
		{"pdb",                1, NULL,               'p'},
		{"prefix",             1, NULL,               'x'},
		{"no-sat-corr",        0, &config_satcorr,     0},
		{"sat-corr",           0, &config_satcorr,     1}, /* Compat */
		{"threshold",          1, NULL,               't'},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{"no-closer-peak",     0, &config_closer,      0},
		{"insane",             0, &config_insane,      1},
		{"image",              1, NULL,               'e'},
		{"basename",           0, &config_basename,    1},
		{"bg-sub",             0, &config_bgsub,       1}, /* Compat */
		{"no-bg-sub",          0, &config_bgsub,       0},

		{"peaks",              1, NULL,                2},
		{"cell-reduction",     1, NULL,                3},
		{"min-gradient",       1, NULL,                4},
		{"record",             1, NULL,                5},
		{"cpus",               1, NULL,                6},
		{"cpugroup",           1, NULL,                7},
		{"cpuoffset",          1, NULL,                8},
		{"hdf5-peaks",         1, NULL,                9},
		{"copy-hdf5-field",    1, NULL,               10},
		{"min-snr",            1, NULL,               11},
		{"min-integration-snr",1, NULL,               12},
		{"tolerance",          1, NULL,               13},
		{"int-radius",         1, NULL,               14},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:wp:j:x:g:t:o:b:e:",
	                        longopts, NULL)) != -1)
	{
		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			filename = strdup(optarg);
			break;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 'z' :
			indm_str = strdup(optarg);
			break;

			case 'p' :
			pdb = strdup(optarg);
			break;

			case 'x' :
			prefix = strdup(optarg);
			break;

			case 'j' :
			nProcesses = atoi(optarg);
			break;

			case 'g' :
			geometry = strdup(optarg);
			break;

			case 't' :
			threshold = strtof(optarg, NULL);
			break;

			case 'b' :
			beam = get_beam_parameters(optarg);
			if ( beam == NULL ) {
				ERROR("Failed to load beam parameters"
				      " from '%s'\n", optarg);
				return 1;
			}
			break;

			case 'e' :
			element = strdup(optarg);
			break;

			case 2 :
			speaks = strdup(optarg);
			break;

			case 3 :
			scellr = strdup(optarg);
			break;

			case 4 :
			min_gradient = strtof(optarg, NULL);
			break;

			case 5 :
			stream_flags = parse_stream_flags(optarg);
			if ( stream_flags < 0 ) return 1;
			break;

			case 6 :
			cpu_num = strtol(optarg, &endptr, 10);
			if ( !( (optarg[0] != '\0') && (endptr[0] == '\0') ) ) {
				ERROR("Invalid number of CPUs ('%s')\n",
				      optarg);
				return 1;
			}
			break;

			case 7 :
			cpu_groupsize = strtol(optarg, &endptr, 10);
			if ( !( (optarg[0] != '\0') && (endptr[0] == '\0') ) ) {
				ERROR("Invalid CPU group size ('%s')\n",
				      optarg);
				return 1;
			}
			if ( cpu_groupsize < 1 ) {
				ERROR("CPU group size cannot be"
				      " less than 1.\n");
				return 1;
			}
			break;

			case 8 :
			cpu_offset = strtol(optarg, &endptr, 10);
			if ( !( (optarg[0] != '\0') && (endptr[0] == '\0') ) ) {
				ERROR("Invalid CPU offset ('%s')\n",
				      optarg);
				return 1;
			}
			if ( cpu_offset < 0 ) {
				ERROR("CPU offset must be positive.\n");
				return 1;
			}
			break;

			case 9 :
			hdf5_peak_path = strdup(optarg);
			break;

			case 10 :
			add_copy_hdf5_field(copyme, optarg);
			break;

			case 11 :
			min_snr = strtof(optarg, NULL);
			break;

			case 12 :
			min_int_snr = strtof(optarg, NULL);
			break;

			case 13 :
			toler = strdup(optarg);
			break;

			case 14 :
			intrad = strdup(optarg);
			break;

			case 0 :
			break;

			default :
			return 1;

		}

	}

	if ( (cpu_num > 0) && (cpu_num % cpu_groupsize != 0) ) {
		ERROR("Number of CPUs must be divisible by"
		      " the CPU group size.\n");
		return 1;
	}

	if ( filename == NULL ) {
		filename = strdup("-");
	}
	if ( strcmp(filename, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(filename, "r");
	}
	if ( fh == NULL ) {
		ERROR("Failed to open input file '%s'\n", filename);
		return 1;
	}
	free(filename);

	if ( outfile == NULL ) {
		outfile = strdup("-");
	}
	if ( strcmp(outfile, "-") == 0 ) {
		ofh = stdout;
	} else {
		ofh = fopen(outfile, "w");
	}
//printf("***************%s\n",outfile);
	if ( ofh == NULL ) {
		ERROR("Failed to open output file '%s'\n", outfile);
		return 1;
	}

	if ( hdf5_peak_path == NULL ) {
		hdf5_peak_path = strdup("/processing/hitfinder/peakinfo");
	}

	if ( speaks == NULL ) {
		speaks = strdup("zaef");
		STATUS("You didn't specify a peak detection method.\n");
		STATUS("I'm using 'zaef' for you.\n");
	}
	if ( strcmp(speaks, "zaef") == 0 ) {
		peaks = PEAK_ZAEF;
	} else if ( strcmp(speaks, "hdf5") == 0 ) {
		peaks = PEAK_HDF5;
	} else {
		ERROR("Unrecognised peak detection method '%s'\n", speaks);
		return 1;
	}
	free(speaks);

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
	}

	if ( nProcesses == 0 ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	if ( (indm_str == NULL) ||
	     ((indm_str != NULL) && (strcmp(indm_str, "none") == 0)) ) {
		STATUS("Not indexing anything.\n");
		indexer_needs_cell = 0;
		reduction_needs_cell = 0;
		indm = NULL;
		cellr = CELLR_NONE;
	} else {
		if ( indm_str == NULL ) {
			STATUS("You didn't specify an indexing method, so I "
			       " won't try to index anything.\n"
			       "If that isn't what you wanted, re-run with"
			       " --indexing=<method>.\n");
			indm = NULL;
			indexer_needs_cell = 0;
		} else {
			indm = build_indexer_list(indm_str,
			                          &indexer_needs_cell);
			if ( indm == NULL ) {
				ERROR("Invalid indexer list '%s'\n", indm_str);
				return 1;
			}
			free(indm_str);
		}

		reduction_needs_cell = 0;
		if ( scellr == NULL ) {
			STATUS("You didn't specify a cell reduction method, so"
			       " I'm going to use 'reduce'.\n");
			cellr = CELLR_REDUCE;
			reduction_needs_cell = 1;
		} else {
			int err;
			cellr = parse_cell_reduction(scellr, &err,
			                             &reduction_needs_cell);
			if ( err ) {
				ERROR("Unrecognised cell reduction '%s'\n",
			              scellr);
				return 1;
			}
			free(scellr);
		}
	}

	/* No indexing -> no reduction */
	if ( indm == NULL ) reduction_needs_cell = 0;

	if ( toler != NULL ) {
		int ttt;
		ttt = sscanf(toler, "%f,%f,%f,%f",
		             &tols[0], &tols[1], &tols[2], &tols[3] );
		if ( ttt != 4 ) {
			ERROR("Invalid parameters for '--tolerance'\n");
			return 1;
		}
	}

	if ( intrad != NULL ) {
		int r;
		r = sscanf(intrad, "%f,%f,%f", &ir_inn, &ir_mid, &ir_out);
		if ( r != 3 ) {
			ERROR("Invalid parameters for '--int-radius'\n");
			return 1;
		}
	} else {
		STATUS("WARNING: You did not specify --int-radius.\n");
		STATUS("WARNING: I will use the default values, which are"
		       " probably not appropriate for your patterns.\n");
	}

	if ( geometry == NULL ) {
		ERROR("You need to specify a geometry file with --geometry\n");
		return 1;
	}

	det = get_detector_geometry(geometry);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);

	if ( reduction_needs_cell || indexer_needs_cell ) {
		cell = load_cell_from_pdb(pdb);
		if ( cell == NULL ) {
			ERROR("Couldn't read unit cell (from %s)\n", pdb);
			return 1;
		}
	} else {
		STATUS("No cell needed for these choices of indexing"
		       " and reduction.\n");
		cell = NULL;
	}
	free(pdb);
	
	write_stream_header(ofh, argc, argv);

	if ( beam != NULL ) {
		nominal_photon_energy = beam->photon_energy;
	} else {
		STATUS("No beam parameters file was given, so I'm taking the"
		       " nominal photon energy to be 2 keV.\n");
		nominal_photon_energy = 2000.0;
	}

	/* Get first filename and use it to set up the indexing */
	prepare_line = malloc(1024*sizeof(char));
	rval = fgets(prepare_line, 1023, fh);
	if ( rval == NULL ) {
		ERROR("Failed to get filename to prepare indexing.\n");
		return 1;
	}
	chomp(prepare_line);
	if ( config_basename ) {
		char *tmp;
		tmp = safe_basename(prepare_line);
		free(prepare_line);
		prepare_line = tmp;
	}
	snprintf(prepare_filename, 1023, "%s%s", prefix, prepare_line);
	qargs.use_this_one_instead = prepare_line;

	/* Prepare the indexer */
	if ( indm != NULL ) {
		ipriv = prepare_indexing(indm, cell, prepare_filename, det,
		                         nominal_photon_energy);
		if ( ipriv == NULL ) {
			ERROR("Failed to prepare indexing.\n");
			return 1;
		}
	} else {
		ipriv = NULL;
	}

	gsl_set_error_handler_off();

	qargs.static_args.cell = cell;
	qargs.static_args.config_cmfilter = config_cmfilter;
	qargs.static_args.config_noisefilter = config_noisefilter;
	qargs.static_args.config_verbose = config_verbose;
	qargs.static_args.config_polar = config_polar;
	qargs.static_args.config_satcorr = config_satcorr;
	qargs.static_args.config_closer = config_closer;
	qargs.static_args.config_insane = config_insane;
	qargs.static_args.config_bgsub = config_bgsub;
	qargs.static_args.cellr = cellr;
	qargs.static_args.tols[0] = tols[0];
	qargs.static_args.tols[1] = tols[1];
	qargs.static_args.tols[2] = tols[2];
	qargs.static_args.tols[3] = tols[3];
	qargs.static_args.threshold = threshold;
	qargs.static_args.min_gradient = min_gradient;
	qargs.static_args.min_snr = min_snr;
	qargs.static_args.min_int_snr = min_int_snr;
	qargs.static_args.det = det;
	qargs.static_args.indm = indm;
	qargs.static_args.ipriv = ipriv;
	qargs.static_args.peaks = peaks;
	qargs.static_args.output_mutex = &output_mutex;
	qargs.static_args.ofh = ofh;
	qargs.static_args.beam = beam;
	qargs.static_args.element = element;
	qargs.static_args.stream_flags = stream_flags;
	qargs.static_args.hdf5_peak_path = hdf5_peak_path;
	qargs.static_args.copyme = copyme;
	qargs.static_args.nProcesses = nProcesses;
	qargs.static_args.ir_inn = ir_inn;
	qargs.static_args.ir_mid = ir_mid;
	qargs.static_args.ir_out = ir_out;

	qargs.fh = fh;
	qargs.prefix = prefix;
	qargs.config_basename = config_basename;
	qargs.n_indexable = 0;
	qargs.n_processed = 0;
	qargs.n_indexable_last_stats = 0;
	qargs.n_processed_last_stats = 0;
	qargs.n_indexableTotal = 0;
	qargs.n_processedTotal = 0;
	qargs.n_indexable_last_statsTotal = 0;
	qargs.n_processed_last_statsTotal = 0;
	qargs.t_last_stats = get_time();
	qargs.updateReader = nProcesses-1; // last process reads updates	

	//////////// Read .lst file ///////////////
	int i,j;
	double t1, t2;
	rewind(fh); // make sure count from start
    n_images = count_patterns(fh);
    rewind(fh);
	//  Divide images into nProcesses
    int nPerBatch = (int)ceil((float)n_images/nProcesses);
    // Malloc 3D array
    char ***array;
    int nChars = LINE_LENGTH;
    array = malloc(nProcesses * sizeof(char **));
    if (array == NULL){
        printf("Error\n");
        return -1;
    }
    for(i=0;i<nProcesses;i++){
        array[i] = malloc(nPerBatch * sizeof(char *));
        if(array[i] == NULL){
            printf("Error\n");
            return -1;
        }
        for(j=0;j<nPerBatch;j++){
            array[i][j] = malloc(nChars * sizeof(char));
            if(array[i][j] == NULL){
                printf("Error\n");
                return -1;
            }
        }
    }
    // Read from file
    for (i=0; i<nProcesses; i++) {
        for (j=0; j<nPerBatch; j++) {
            if(fgets(&array[i][j][0],(nChars * sizeof(char))-1,fh) == NULL){
				array[i][j] = NULL;
			}
        }
    }
    fclose(fh);
    printf("\n"); // flush buffer

	// Clear output file content
	char *myOutfilename = NULL;
	chomp(prefix);
	chomp(outfile);
	myOutfilename = malloc(strlen(outfile)+1);
	snprintf(myOutfilename,LINE_LENGTH-1,"%s",outfile); //maximum print length
	FILE *tfh;	
	tfh = fopen(myOutfilename,"a+");
	if (tfh == NULL) { printf("Error\n"); }
	fclose(tfh);
	qargs.static_args.outfile = outfile;	

	//////////// PIPING ///////////////
	char bufferW[BUFSIZE], bufferR[BUFSIZE]; // write and read buffers
	int fd_pipe[nProcesses][2];	

	memset(bufferW, 0, BUFSIZE);
	memset(bufferR, 0, BUFSIZE);

	for (i=0;i<nProcesses;i++){	// pipe nthread times
		pipe(fd_pipe[i]);
	}

	//////////// FORKING ///////////////
	int power = 10;  // 2^power must be larger than nProcesses
	int pid[power];
	double num = 0;
	int batchNum = 0;
	// Fork 2^power times
	for (i=0;i<power;i++){
		pid[i] = fork();
	}
	// Assign id
	for (i=0;i<power;i++){
		if (pid[i] == 0){ // keep parents and kill off children
			num += pow(2,i);
		}
	}
	// Kill if batchNum too high
	if (num >= nProcesses){
		exit(0); // kill
	}
	batchNum = (int) num;

	// Calculate how many images to process
	qargs.nPerProcess = 0;
	i=batchNum;
    for (j=0; j<nPerBatch; j++) {
		if (array[i][j] != '\0'){
			qargs.nPerProcess++;
		}
    }
	//printf("%d: %d images\n",batchNum,qargs.nPerProcess);

	////////// PLUMBING //////////
	if (batchNum == qargs.updateReader){
		for (i=0;i<nProcesses;i++){
			if (i != batchNum) {
				close(fd_pipe[i][1]); // only open own write pipe
			}
			// read from all processes
		}		
	} else {
		for (i=0;i<nProcesses;i++){
			if (i != batchNum) {
				close(fd_pipe[i][1]); // only open own write pipe
			}
			close(fd_pipe[i][0]); // close all read pipes
		}
	}

	//////////// INDEXING ///////////////
	//printf("PID&PPID: %d %d\n",getpid(),getppid());

	t1 = get_time();

	process_image(array,batchNum,&qargs,fd_pipe);

	///// READ UPDATE /////
	if (batchNum == qargs.updateReader){ // loop until every process is finished
		int pipesOpen;
		int loopFlag = 1;
			while (loopFlag) {
			double seconds;
			seconds = get_time();
			if ( seconds >= qargs.t_last_stats+STATS_EVERY_N_SECONDS ) {
				pipesOpen = readUpdate(&qargs,fd_pipe,1);
				if (pipesOpen == 0) {
					// Close all read pipes
					for (i=0;i<qargs.static_args.nProcesses;i++){
						close(fd_pipe[i][0]); 
					}
					loopFlag = 0; // Must execute this to get out
				}			
			}
		}
	}

	close(fd_pipe[batchNum][1]); // close own write pipe

	t2 = get_time();

	//////////// EXITING ///////////////	

	printf("%d: Fork Time %.2fs\n",batchNum,t2-t1);

	if (batchNum == qargs.updateReader){
		STATUS("There were %i images, of which %i could be indexed.\n",
	        n_images, qargs.n_indexableTotal);
	}

	exit(0);
}
