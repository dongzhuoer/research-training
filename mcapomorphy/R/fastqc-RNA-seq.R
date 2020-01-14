#' @title browse fastqc html output view RStduio Server
#' @section package options embedded in aves:
#'     options(mcapomorphy.wd = '/path/to/wd')
#'     options(mcapomorphy.url = 'http://hostname:port')



#' @title return the url of a file on the workstation
#'
#' @param path string. full absolute path of the file
#'
#' @return string.
#'
#' @examples
#' file_url(fastqc_path('SRR6148275_2'))
#'
#' @keywords internal
file_url <- function(path) {
	paste0(getOption('mcapomorphy.url'), '/file_show?path=', path)
}



#' @title return the full absolute path of the fastqc result
#'
#' @param name string. name of the sequence file, such as `'SRR6148275_2'`
#'
#' @return string.
#'
#' @examples 
#' fastqc_path('SRR6148275_2')
#'
#' @keywords internal
fastqc_path <- function(name) {
	paste0(getOption('mcapomorphy.wd'), '/oases/fastqc/', name, '_fastqc.html')
}



#' @title view fastqc result of original SRA sequence for single ends
#'
#' @param id string. sra accession, such as `'SRR029421'`
#'
#' @return `NULL`
#'
#' @examples
#' \dontrun{
#'     view_sra_single('SRR029421')
#' }
#'
#' @keywords internal
view_sra_single <- function(id) {
	id %>% paste0('output/', .) %>% fastqc_path %>% file_url %>% libzhuoer::browse_url();
}



#' @title view fastqc result of original SRA sequence for paired ends
#'
#' @param id string. sra accession, such as `'SRR6148275'`
#'
#' @return `NULL`
#'
#' @examples
#' \dontrun{
#'     view_sra_paired('SRR6148275')
#' }
#'
#' @keywords internal
view_sra_paired <- function(id) {
	id %>% paste0('output/', .) %>% paste0('_', 1:2) %>% fastqc_path %>% file_url %>% libzhuoer::browse_url();
}



#' @title compare fastqc result of original and trimmed sequence for single ends
#'
#' @param id string. sra accession, such as `'SRR029421'`
#'
#' @return `NULL`
#'
#' @examples
#' \dontrun{
#'     compare_sra_single('SRR029421')
#' }
#'
#' @keywords internal
compare_sra_single <- function(id) {
	id %>% paste0(c('input/', 'output/'), .) %>% fastqc_path %>% file_url %>% libzhuoer::browse_url();
}



#' @title compare fastqc result of original and trimmed sequence for paired ends
#'
#' @param id string. sra accession, such as `'SRR6148275'`
#'
#' @return `NULL`
#'
#' @examples
#' \dontrun{
#'     compare_sra_paired('SRR6148275')
#' }
#'
#' @keywords internal
compare_sra_paired <- function(id) {
	names <- id %>% paste0('_', 1:2) %>% plyr::llply(. %>% paste0(c('input/', 'output/'), .)) %>% unlist;
	names %>% fastqc_path %>% file_url %>% libzhuoer::browse_url(0.5);
}


