function (experiment, annotation, chromosomes, samples, consensus = FALSE, 
	bCount = FALSE, mapQCth = 15, blacklist = NULL, profileWin = 400, 
	fragmentLength = 125, shifts = 1:300, ...) 
{
	if (class(experiment) == "character" || class(experiment) == 
		"data.frame") {
		experiment = dba(sampleSheet = experiment, bCorPlot = FALSE, 
			peakCaller = "bed")
	}
	if (class(experiment) != "DBA") {
		stop("experiment must be either a samplesheet filename or a DBA (DiffBind) object.")
	}
	experiment$config$mapQCth = mapQCth
	meta = data.frame(t(experiment$class))
	if (length(unique(meta$bamRead)) != nrow(meta)) {
		stop("Unable to process. Each bam file must be associated with at most one peakset.")
	}
	if (!missing(samples)) {
		for (i in 1:length(experiment$peaks)) {
			if (nrow(experiment$peaks[[i]]) == 0) {
				experiment = addMatchingSample(experiment, i, 
					meta, samples)
			}
		}
		experiment = dba(experiment, bCorPlot = FALSE)
	}
	if (missing(chromosomes)) {
		chromosomes = 1
	}
	if (is.numeric(chromosomes) && missing(samples)) {
		chrmap = experiment$chrmap
		if (length(chrmap) == 0) {
			warning("No chromosomes specified in peaks, using all.")
			chromosomes = NULL
		}
		else {
			if (max(chromosomes) > length(chrmap)) {
				warning("Specified chromosome number exceeds chromosomes seen in peaks.")
				chromosomes = chromosomes[chromosomes <= length(chrmap)]
			}
			chromosomes = chrmap[chromosomes]
			message("Checking chromosomes:")
			print(chromosomes)
		}
	}
	if (!missing(annotation)) {
		if (!is.null(annotation) && missing(samples)) {
			if (class(annotation) != "list") {
				message("Compiling annotation...")
				annotation = getAnnotation(annotation, AllChr = chromosomes)
			}
			if (annotation$version == "hg19" && missing(blacklist)) {
				blacklist = read.table(file.path(system.file("extdata", 
					package = "ChIPQC"), "blacklist_hg19.bed"), 
					header = TRUE)[, 1:4]
				blacklist = makeGRangesFromDataFrame(blacklist, 
					ignore.strand = TRUE)
				message("Using default blacklist for hg19...")
			}
		}
		else if (class(annotation) == "character") {
			annotation = list(version = annotation)
		}
		else {
			annotation = list(version = "none")
		}
	}
	else {
		annotation = list(version = "none")
	}
	samplelist = NULL
	controlist = NULL
	for (i in 1:nrow(meta)) {
		newrec = NULL
		newrec$peaks = experiment$peaks[[i]]
		if (nrow(newrec$peaks) == 0) {
			newrec$peaks = NULL
		}
		newrec$bam = as.character(meta$bamRead[i])
		samplelist = listadd(samplelist, newrec)
		if (!is.null(meta$bamControl[i])) {
			if (!is.na(meta$bamControl[i])) {
				if (meta$bamControl[i] != "") {
					savenames = names(controlist)
					controlist = listadd(controlist, as.character(meta$bamControl[i]))
					names(controlist) = c(savenames, as.character(meta$Control[i]))
				}
			}
		}
	}
	controlist = controlist[!duplicated(controlist)]
	controls = 0
	for (cfile in controlist) {
		if (!(cfile %in% as.character(meta$bamRead))) {
			newrec = NULL
			newrec$bam = cfile
			newrec$peaks = NULL
			samplelist = listadd(samplelist, newrec)
			controls = controls + 1
		}
	}
	names(samplelist) = unique(c(rownames(meta), names(controlist)))
	addconsensus = FALSE
	addcontrols = FALSE
	setNull <- TRUE
	peaks <- NA
	if (class(consensus) == "GRanges") {
		useConsensus <- TRUE
	}
	else {
		if (class(consensus) == "logical") {
			if (consensus == TRUE) {
				useConsensus <- TRUE
			}
			else {
				useConsensus <- FALSE
			}
		}
		else {
			stop("consensus must be either a logical or a GRanges.")
		}
	}
	if (useConsensus == TRUE) {
		addcontrols <- TRUE
		addconsensus <- TRUE
		if (class(consensus) == "GRanges") {
			peaks <- consensus
			setNull <- FALSE
		}
		else {
			setNull <- TRUE
			peaks <- dba.peakset(experiment, bRetrieve = T, 
				numCols = 3, DataType = DBA_DATA_FRAME)
		}
	}
	else if (bCount) {
		addcontrols <- TRUE
		if (consensus == FALSE) {
			setNull <- TRUE
			peaks <- dba.peakset(experiment, bRetrieve = TRUE, 
				DataType = DBA_DATA_FRAME)[, 1:4]
		}
	}
	if (addcontrols && controls) {
		message("Adding controls...")
		startpos = length(samplelist) - controls
		for (i in 1:controls) {
			metadata = getMeta(meta, names(samplelist)[[startpos + 
				i]])
			experiment = dba.peakset(experiment, peaks = peaks, 
				sampID = names(samplelist)[[startpos + i]], 
				factor = "Control", tissue = metadata$tissue, 
				condition = metadata$condition, treatment = metadata$treatment, 
				replicate = sprintf("c%s", i), bamReads = samplelist[[startpos + 
					i]]$bam)
		}
		meta = data.frame(t(experiment$class))
	}
	if (is.na(peaks) || setNull == TRUE) {
		peaks = NULL
	}
	if (bCount) {
		message("Counting reads in consensus peakset...")
		if (!is.null(peaks)) {
			experiment = dba.count(experiment, peaks = peaks, 
				mapQCth = mapQCth, ...)
			meta = data.frame(t(experiment$class))
		}
		else {
			experiment = dba.count(experiment, mapQCth = mapQCth, 
				...)
			meta = data.frame(t(experiment$class))
		}
		if (nrow(meta) != length(samplelist)) {
			warning("Samples and peaksets out of sync")
		}
	}
	if (nrow(experiment$merged) > 0) {
		peaks = dba.peakset(experiment, bRetrieve = TRUE)
	}
	else peaks = NULL
	if (addcontrols && addconsensus) {
		for (i in 1:length(samplelist)) samplelist[[i]]$peaks = peaks
	}
	else if (addcontrols) {
		for (i in 1:length(samplelist)) {
			if (is.null(samplelist[[i]]$peaks)) {
				samplelist[[i]]$peaks = peaks
			}
		}
	}
	if (missing(samples)) {
		message(sprintf("Computing metrics for %d samples...", 
			length(samplelist)))
		samples = bplapply(samplelist, doChIPQCsample, experiment, 
			chromosomes, annotation, mapQCth, blacklist, profileWin, 
			fragmentLength, shifts)
		names(samples) = names(samplelist)
		errors = FALSE
		for (i in 1:length(samples)) {
			sample = samples[[i]]
			if (class(sample) != "ChIPQCsample") {
				message("Error in sample ", names(samples)[i], 
					": ", sample[1])
				errors = TRUE
			}
		}
		if (errors) {
			stop("Errors in generating sample data.")
		}
	}
	else {
		if (length(samples) < length(samplelist)) {
			stop("Not enough samples!")
		}
		else if (length(samples) > length(samplelist)) {
			warning("Ignoring some samples without metadata.")
		}
		newsamps = NULL
		for (sampname in names(samplelist)) {
			whichsamp = which(names(samples) %in% sampname)
			if (length(whichsamp) == 0) {
				stop(sprintf("Sample %d missing from sample list", 
					sampname))
			}
			else if (length(whichsamp) > 1) {
				stop(sprintf("Sample %d appears more than once in sample list", 
					sampname))
			}
			newsamps = listadd(newsamps, samples[[whichsamp]])
			names(newsamps)[length(newsamps)] = sampname
		}
		samples = newsamps
	}
	res = new("ChIPQCexperiment", Samples = samples, DBA = experiment, 
		annotation = annotation)
	config = as.list(res@DBA$config)
	config$fragmentSize = fragmentlength(res)
	config$mapQCth = mapQCth
	res@DBA$config = config
	return(res)
}
