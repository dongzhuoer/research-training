#' @title wrapper to run AUGUSTUS
#'
#'
#' @param species string. an identifier for the species model to use, such as `'chicken'`
#' @param input string. path to the input file, which contains the gene sequence(s) in fasta format.
#' @param output string. path to the output file. if `NULL` output is writen into standard output
#' @param args character. additional parameters passed on to AUGUSTUS
#' @param ... other parameters passed on to [system2()], from `stdout` to the end
#'
#' @return result of [system2()] or `NULL` if `input` doesn't exist
#'
#' @examples
#' \dontrun{
#'     augustus('human', 'software/augustus/examples/example.fa');
#'     augustus('human', 'software/augustus/examples/hsackI10.fa.gz')
#'
#'     augustus('human', 'software/augustus/examples/example.fa', 'test.gff');
#' }
#' 
#'
#' @export
augustus <- function(species, input, output = NULL, args = '', ...) {
	if (!file.exists(input)) return(NULL);

	species %<>% paste0('--species=', .);
	output %<>% {if (is.null(.)) '' else paste0('--outfile=', .)};
	args %<>% c(species, input, output, .);

	system2('augustus', args, ...);
}






#' go to http://www.orthodb.org/, run `$('#full-tree').fancytree('getTree').visit(function(node){node.setExpanded()})` in console, save page all
#' search 'fancytree-has-children' for internal node, 'fancytree-checkbox' for all node
#'


# get_children <- fucntion(node){
# 	has_children <- is.list(node) && 'children' %in% names(node);
# 	if (has_children) node$children else list();
# }
#
# orthodb_all_node <- function() {
# 	pkg_file('extdata/orthodb-tree.json.gz') %>% jsonlite::read_json() %>%
# 		{lapply(.$data, libzhuoer::nested_element, get_children)} %>% unlist(F)
# }
#
# orthodb_internal_node <- function() {
# 	orthodb_all_node() %>% {.[is_internal(.)]}
# }
#
# is_internal_node <- . %>% sapply(. %>% names %>% is.element('children', .));
#
# ls_orthodb_clade <- function() {
# 	orthodb_internal_node() %>% sapply(. %>% {.$name})
# }
#
# get_species <- function(clade_name) {
# 	orthodb_internal_node()[[which(ls_orthodb_clade() == clade_name)]] %>%
# 		libzhuoer::nested_element(get_node) %>% {.[!is_internal(.)]} %>%
# 		lapply(tibble::as_tibble) %>% dplyr::bind_rows()
# }
# get_species('Aves')
#
# readr::read_lines('http://www.orthodb.org/tree') %>% readr::write_lines(gzfile('inst/extdata/orthodb-tree.json.gz'))
# jsonlite::read_json('data-raw/orthodb-tree.json.gz')
#
#
# orthodb_tree <- jsonlite::read_json('data-raw/orthodb-tree.json.gz')
#
#
# all_node <-
#
#
# is_internal(all_node) %>% sum
#
# all_node %>% {.[is_internal(.)]} %>% sapply(. %>% {.$name})
#
#
#
# internal_node <- all_node[is_internal]
# tip_node <- all_node[!is_internal]
#
# sapply(internal_node, . %>% {.$name})











