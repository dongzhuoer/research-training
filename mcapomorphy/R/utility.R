#' @title get files in the package
#' @export
pkg_file <- function(..., mustWork = TRUE) {
    system.file(..., package = "mcapomorphy", mustWork = mustWork)
}



#' @title replace extension name
#'
#' @description extension name can be any number of `\\w`
#'
#' @examples
#' replace_ext('a.txt', '.md')
#'
#' @section to do: unittest
#'
#' @keywords internal
replace_ext <- function(filename, new_ext = '') {
	stringr::str_replace(filename, '\\.\\w*$', new_ext);
}




#' @title get the outermost list which contains more than one element
#'
#' @description `unnest_list(list(list(list(1))))` gives `list(1)`
#'
#' @examples
#' unnest_list(1)
#' unnest_list(list(1))
#' unnest_list(list(list(list(1))))
#' unnest_list(list(list(list(1), list(2))))
#' unnest_list(list(list(list(1)), list(2)))
#'
#' @export
unnest_list <- function(x) {
	if (!is.list(x)) return(x);

	while (is.list(x[[1]]) && length(x) == 1L) {
		x = x[[1]]
	}
	x
}








