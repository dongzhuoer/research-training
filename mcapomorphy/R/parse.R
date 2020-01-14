#' @title parse NCBI genome webpage
#'
#' @param path string. path to the HTML file.
#'
#' @return [tibble::tibble()]
#'
#'   1. species
#'
#'   1. number
#'
#'   1. dna
#'
#'   1. rna
#'
#'   1. protein
#'
#' @examples
#' parse_genome(system.file('extdata', 'genome55342.html', package = 'mcapomorphy'))
#' parse_genome(system.file('extdata', 'genome2793.html', package = 'mcapomorphy'))
#'
#' @export
parse_genome <- function(path) {
	html <- xml2::read_html(path);

	species <- html %>% rvest::html_nodes('span.GenomeTitle') %>%
		rvest::html_text() %>% stringr::str_extract('\\w+ \\w+');

	number <- html %>% rvest::html_nodes('.refgenome_sensor b') %>%
		rvest::html_text() %>% stringr::str_subset('\\d+ genomes') %>%
		{if (length(.) > 0) stringr::str_extract(., '\\d+') %>% as.integer else 1L};

	get_url <- function(pattern) {
		html %>% rvest::html_nodes('span.shifted a') %>% rvest::html_attr('href') %>%
			stringr::str_subset(pattern) %>% {if (length(.) > 0) . else ''}
	}

	tibble::data_frame(species = species, number = number,
					   dna = get_url('genomic.fna'), rna = get_url('rna.fna'),
					   protein = get_url('protein.faa'));
}



#' @title parse assembly metadata in xml format
#'
#' @param path `x` of [xml2::read_xml()]. assembly metadata in XML format
#'
#' @return [tibble::tibble()]
#'
#'   1. species
#'
#'   1. id
#'
#'   1. similarity logical scalar.
#'
#'   1. partial logical scalar.
#'
#'   1. gca_dir
#'
#'   1. gcf_dir
#'
#' @section to do: unittest
#'
#' @examples
#' NULL
#'
#' @export
parse_assembly <- function(path) {
	as_string <- . %>% {if (rlang::is_empty(.)) '' else .[[1]]};

	assembly <- xml2::read_xml(path) %>% xml2::as_list() %>% `[[`(1);
	species = assembly$SpeciesName[[1]];
	id = assembly$AssemblyAccession[[1]];
	similarity = assembly$Synonym$Similarity %>%
		{if (rlang::is_empty(.)) NA else .[[1]]};
	partial = as.logical(assembly$PartialGenomeRepresentation[[1]]);
	gca_dir = as_string(assembly$FtpPath_GenBank);
	gcf_dir = as_string(assembly$FtpPath_RefSeq);

	tibble::tibble(species, id, similarity, partial, gca_dir, gcf_dir);
}



#' @title parse avibase species webpage
#'
#' @details you have to download a Chinese webpage: note `lang=ZH` in `https://avibase.bsc-eoc.org/species.jsp?lang=ZH&avibaseid=`
#'
#' @param path string. path to the HTML file.
#'
#' @return data_frame
#'
#' 1. order
#'
#' 1. family
#'
#' 1. genus
#'
#' 1. species
#'
#' 1. id
#'
#' @examples
#' parse_avibase(system.file('extdata', 'Avibase-01761045A529DC4A.html', package = 'mcapomorphy'))
#' parse_avibase(system.file('extdata', 'Avibase-1FDDABDB0D4421F9.html', package = 'mcapomorphy'))
#'
#' @export
parse_avibase <- function(path) {
	avibase <- path %>% xml2::read_html();

	species <- avibase %>% rvest::html_node('td.AVBContainerText p i') %>% rvest::html_text();
	v <- avibase %>% rvest::html_nodes('td.AVBContainerText script + p') %>%
		rvest::html_text() %>% {stringr::str_extract_all(., '[a-zA-Z]+')[[1]]};

	tibble::data_frame(order = v[1], family = v[2], genus = v[3],
					   species = species, id = basename(path));
}
