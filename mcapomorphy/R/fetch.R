
#' @title download NCBI genome webpage for a vector of given id
#'
#' @param assessions character.
#'
#' @return `NULL`
fetch_genome <- function(id) {
	fetch(id, 'https://www.ncbi.nlm.nih.gov/genome/', 'data-raw/genome-html/');
}



#" for debug `fetch()`
# srcs = c('45123', '66440', '61239', '58928', '53567', '50144', '46404', '43765', '40104');
# prefix_url = 'https://www.ncbi.nlm.nih.gov/genome/';
# dest_dir = 'data-raw/genome-html/';

#' @title download many webpages parallelly
#'
#' @details only retry for at most `getOption('mcapomorphy.max_iterate', 64L)`
#'   times per execute, you may run it again if you still have problems
#'
#' @param srcs character. would become the names of downloaded files.
#' @param prefix_url string. such as `'https://en.wikipedia.org/wiki/'`,
#'   `paste0(prefix_url, query)` gives the url, see `mutate` for `query`
#' @param dest_dir string. such as 'data-raw/wikipedia/'
#'
#' @return `NULL`
#' @export
fetch <- function(srcs, prefix_url, dest_dir) {
	dest_dir %>% {if (!dir.exists(.)) dir.create(.)};
	remove_empty <- . %>% {if (file.size(.) == 0 ) file.remove(.)};
	i <- 0L;
    max_iterate <- getOption('mcapomorphy.max_iterate', 256L);

	while (TRUE) {
		dir(dest_dir, full.names = T) %>% plyr::l_ply(remove_empty);

	    done <- dest_dir %>% dir;
	    done = done[file.size(paste0(dest_dir, done)) > 10];
	    todo = setdiff(srcs, done);

	    i = i + 1;
	    if (i >= max_iterate) {
	    	warning('reached max iterate, stoping');
	    	break;
	    }

	    if (length(todo) == 0) break;

	    parallel::mclapply(
	    	todo,
	    	function(src) {
	    		query <- src;

				url  <- paste0(prefix_url, query);
				dest <- paste0(dest_dir, src);

				message(url, dest)

                tryCatch(download.file(url, dest), error = function(e){});
	    	}
	    );
	}

	NULL;
}
