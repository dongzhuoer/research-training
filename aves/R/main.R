

#' @title command for running trimmomatic
#'
#' @param id string. SRA accession
#' @param adapter string. `'TruSeq3'` or `'TruSeq2'`
#' @param extra 
#'
#' @return string. 
#'
#' @section Adapter: according to the "Overrepresented Sequences" report of FASTQC
#' 
#' 1. "Illumina Single End" or "Illumina Paired End" : TruSeq2
#' 1. "TruSeq Universal Adapter" or "TruSeq Adapter, Index ..." : TruSeq3 
#' 
#' @examples
#' 
#' @export
trimmomatic_command <- function(id, adapter = 'TruSeq3', extra = '') {
	mode <- '';
	input <- '';
	output <- '';
	if (paste0('raw/fastq/', id, '_', 1:2, '.fastq.gz') %>% file.exists %>% all) {
		mode = 'PE'
		input = paste0('raw/fastq/',    id, '_', 1:2, '.fastq.gz', collapse = ' ');
		output = paste0('oases/input/', id, '_', 1:2, '.fastq.gz /dev/null', collapse = ' ');
	} else if (paste0('raw/fastq/', id,'.fastq.gz') %>% file.exists) {
		mode   = 'SE';
		input  = paste0('raw/fastq/',   id,'.fastq.gz');
		output = paste0('oases/input/', id,'.fastq.gz');
	} else {
		warning(id, ' not found!');
		return('');
	}
    clip <- paste0('ILLUMINACLIP:software/Trimmomatic-0.36/adapters/', adapter, '-', mode, '.fa:2:30:10');
	
	paste('trimmomatic', mode, input, output, clip, 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:3:15', 'MINLEN:36');
}





# oases_command <- function() {
# 	[ -d oases/@ ] || mkdir oases/@; nohup time oases_pipeline.py -m 21 -M 45 -s 4 -o oases/@/ -d '-fastq.gz oases/input/@-trim.fastq.gz' &>　oases/log/@ &
# [ -d oases/@ ] || mkdir oases/@; nohup time oases_pipeline.py -m 21 -M 23      -o oases/@/ -d '-fastq.gz -shortPaired -separate oases/input/@_1-trim-paired.fastq.gz oases/input/@_2-trim-paired.fastq.gz' &>　oases/log/@ &
# 
# }



#' @title return ref species for a species
#'
#' @param name character. see [extract_genus]
#'
#' @return character. ref_species paratmeter for HaMStR
#' @export
#'
#' @examples {
#'     aves::ref_species(aves::transcript$name);
#'     aves::ref_species(aves::outgroup$name);
#' }
ref_species <- function(name) {
	genus <- extract_genus(name) %>% 
		{ifelse(. %in% aves::taxonomy$genus, ., 'Gallus')};
	#" let outgroup use Gallus gallus
	order <- aves::taxonomy %>% {.$order[match(genus, .$genus)]};
	
	impl <- function(order) {
		if (order %in% aves::stratum$Telluraves) return('Taeniopygia_guttata')
		if (order %in% aves::stratum$`Sup Aequornithia`) return('Nipponia_nippon')
		if (order %in% aves::stratum$`Other Neoaves`) return('Taeniopygia_guttata,Nipponia_nippon -relaxed')
		if (order %in% aves::stratum$Galloanseres) return('Gallus_gallus')
		if (order %in% aves::stratum$Palaeognathae) return('Gallus_gallus')
		stop('don\'t know how to specify reference species for ', order);
	}
	
	sapply(order, impl);
}




#' @title generate hamstr command
#'
#' @param taxon character. taxon name, see [extract_genus]
#' @param dna logical scalar. Whether the input sequence is EST or protein
#' @param dir string. folder where the input file lies and output file will be
#'   put
#' @param hmmer_top integer scalar. passed on to `hit_limit` option. specify
#'   `NULL` to disable that option.
#' @param hmm_set character.  passed on to `hmmset` option.
#' @param threads integer scalar. how many tasks to run in parallel
#'
#' @return character
#' @section to do: remove `mod` `mod.tc` if older
#' @export
#'
#' @examples {
#'     aves::hamstr_command('Arborophila_rufipectus_SRA_1')
#'     aves::hamstr_command('Gallus_Gallus', F, 'hamstr-ext', NULL, 'Gallus_Gallus')
#' }
hamstr_command <- function(taxon, hmmpath = '/path/to/data/aves/hamstr', hmmer_top = 10L, dir = 'hamstr', 
 hmm_set = 'aves', threads = getOption("mc.cores")) {
	if (length(taxon) == 0L) return('')
	
    input_type <- ifelse(str_detect(taxon, 'genome|outgroup|protein'), ' -protein', ' -est');
    input_postfix <- ifelse(str_detect(taxon, 'genome|outgroup|protein'), '.faa', '.fna');
    cpu <- if (is.null(threads)) '' else paste0(' -cpu ', threads);
    hit_limit <- if (is.null(hmmer_top)) '' else paste0(' -hit_limit ', hmmer_top);
    
    paste0('echo; echo ', taxon, '; time hamstr', input_type, ' -outpath=', dir, ' -sequence_file=', dir, '/', taxon, input_postfix, ' -taxon=', taxon, ' -blastpath /path/to/data/aves/hamstr-data/blast_dir', ' -hmmpath=', hmmpath, ' -hmmset=', hmm_set, ' -refspec=', ref_species(taxon), ' -representative -reuse -central -force', cpu, hit_limit, ' &> ', dir, '/', taxon, '.log; sleep 60;')
}








