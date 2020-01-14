#' @title convert fasta to nexus data block.
#'
#' @description internal function for [write_paup()]
#'
#' @details the return value will be printed when call at console, if the file
#'   is big, it will take a while, and you can't do anything. Consider
#'   `fasta_to_nexus() -> x` to avoid automatically printing.
#'
#'   Alternatively, you can use seaview But it may take several seconds if each
#'   sequence lies on only one line in the original file, because seaview always
#'   force the line width less than 80.
#'
#' @param fasta tibble. see **Value** of [biozhuoer::read_fasta()]. Sequence
#'   must be aligned (in the same length). `'U'` would be replaced by
#'   `'X'` since paup doesn't recognized it
#' @param datatype string. what kind of data is you fasta. (`'dna'`, `'rna'`,
#'   `'protein'`, or `'standard'`)
#' @param missing string. symbol used for missing data
#' @param gap string. symbol used for gap
#'
#' @return character. A Nexus data block. Note that begining '#nexus' isn't
#'   included.
#'
#' @examples
#' fasta <- biozhuoer::read_fasta(system.file('extdata', 'EOG090F05Z3.fasta', package = 'mcapomorphy'))
#' fasta_to_nexus(fasta, 'protein')
#' 
#' @section to do: find a OG containing 'U' to testthat
#'
#' @keywords internal
fasta_to_nexus <- function(fasta, datatype, missing='?', gap='-') {
	fasta$name %<>% {paste0(format(.), '    ')};   #add extra spaces to align name
    fasta$seq  %<>% stringr::str_replace_all('U', 'X');

	n <- nrow(fasta);
	if (n == 0L) stop('empty fasta');

	len <- nchar(fasta$seq) %>% unique
	if (length(len) > 1L) stop('unaligned fasta');

	nexus <- character(6 + n);
	nexus[1] = 'begin data;';
	nexus[2] = paste0('dimensions ntax=', n, ' nchar=', len, ';');
	nexus[3] = paste0('format datatype=', datatype, ' missing=', missing, ' gap=', gap, ';');
	nexus[4] = 'matrix';
	nexus[5:(4 + n)] = paste0(fasta$name, fasta$seq);
	nexus[5 + n] = ';';
	nexus[6 + n] = 'end;';

	nexus
}





#' @title perpare nexus file to run paup
#'
#' @details for `fasta_file`, `newick_file` and `outgroup` only
#'   \code{[[:alnum:]_]} are allowed in taxon name (sequence header)
#'
#' @param fasta_file string. sequence file in FASTA format. see
#'   [fasta_to_nexus()] for requirement of the content
#' @param newick_file string. tree file in Newick format.
#' @param outgroup character. name of outgroups, `' '` would be replaced by
#'   `'_'`
#' @param output_file string. path to the output NEXUS file.
#'
#' @return `NULL`
#'
#' @examples
#' temp_file <- tempfile()
#' write_paup(
#'     system.file('extdata', 'EOG090F05Z3.fasta', package = 'mcapomorphy'), 
#'     system.file('extdata', 'omics.tre', package = 'mcapomorphy'), 
#'     aves::outgroup$species, temp_file
#' )
#' \dontrun{
#'     system2('paup4a159', 'data-raw/EOG090F05Z3.nexus -n -u', T)
#' }
#' 
#' @export
write_paup <- function(fasta_file, newick_file, outgroup, output_file) {
	sequence <- biozhuoer::read_fasta(fasta_file) %>% fasta_to_nexus('protein');

	tree <- readr::read_lines(newick_file) %>% c('begin trees;', 'Tree my_tree=', ., 'end;');

	outgroup %<>% stringr::str_replace_all(' ', '_') %>% paste(., collapse = ' ');
	command <- c('begin paup;',
	  paste('outgroup', outgroup, '/only'),
	  'pset opt=acctran',
	  'describetrees /plot=both labelnode=yes apolist=yes',
	  'quit',
	  'end'
	) %>% paste0(';')

	c('#nexus',  sequence, tree, command) %>% readr::write_lines(output_file);

	return(NULL);
}















