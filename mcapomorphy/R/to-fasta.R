
#' @title prepare hamstr input file from augustus gff
#'
#' @details convert gff file produced by augustus to fasta file storing protein
#'
#' @param input. string. path to input file
#' @param output string. path to output file
#'
#' @return `NULL`
#' @export
#'
hamstr_from_augustus <- function(input, output) {
	if (!file.exists(input)) return(NULL);

    gff <- readr::read_lines(input);

	seq <- gff %>% stringr::str_subset('^#') %>%
    	stringr::str_replace('^# ', '') %>%
		paste0(collapse = '') %>%
		{stringr::str_extract_all(., '(?<=protein sequence = \\[)[^\\[\\]]+(?=\\])')[[1]]}
	name <- seq_along(seq);
	fasta <- tibble::tibble(name, seq);

    biozhuoer::write_fasta(fasta, output);

	return(NULL);
}




#' @title prepare hamstr input file from ncbi protein
#'
#' @param input. string. path to input file
#' @param output string. path to output file
#'
#' @return `NULL`
#' @export
#'
hamstr_from_ncbi_protein <- function(input, output) {
	if (!file.exists(input)) return(NULL);

	biozhuoer::read_fasta(input) %>%
        mutate(name = str_replace(name, '\\.[\\w\\W]+$', '')) %>%
        biozhuoer::write_fasta(output);

	return(NULL);
}


#' @title prepare hamstr input file from ncbi TSA
#'
#' @param input. string. path to input file
#' @param output string. path to output file
#'
#' @return `NULL`
#' @export
#'
hamstr_from_TSA <- function(input, output) {
	if (!file.exists(input)) return(NULL);

	fasta_raw <- readr::read_lines(input);

    header_line <- stringr::str_detect(fasta_raw, '^>') %>% which;
    new_header <- fasta_raw[header_line] %>%
    	stringr::str_extract('>[\\w\\.]+') %>% stringr::str_replace('\\.', '_')
    fasta_raw[header_line] = new_header;

    readr::write_lines(fasta_raw, output);

	return(NULL);
}





#' @title prepare hamstr input file from oases transcript
#'
#' @param input. string. path to input file
#' @param output string. path to output file
#'
#' @return `NULL`
#' @export
#'
hamstr_from_oases <- function(input, output) {
	if (!file.exists(input)) return(NULL);

    fasta_raw <- readr::read_lines(source);
    fasta_raw %<>% stringr::str_replace('/[^/]+$', '') %>%
    	stringr::str_replace_all('(Locus|Transcript)_', '');
    readr::write_lines(fasta_raw, dest);

    return(NULL);
}





#' @title convert hamstr output file to FASTA format
#'
#' @details if a OG exists multiple times, only the first one would be reserved
#'
#' @param input. character. path to input file, can be more than one.
#'   non-existent ones are ignored.
#' @param output string. path to output file
#'
#' @return `NULL`
#' @export
#'
hamstr_out_to_fasta <- function(input, output) {
	content <- input %>% {.[file.exists(.)]} %>% lapply(readr::read_lines) %>% unlist;

	fasta <- tibble::data_frame(
		name = stringr::str_extract(content, '^\\w+'),
		seq  = stringr::str_replace(content, '^[\\w\\W]+\\|', '')
	) %>% plyr::ddply('name', . %>% dplyr::slice(1));

    biozhuoer::write_fasta(fasta, output);
}







