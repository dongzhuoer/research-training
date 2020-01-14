#' @title extract the first word (only consists of letters)
#'
#' @description used to extract genus
#'
#' @param name character. can be scientific name or `name` of aves's data.frame
#'
#' @return character
#'
#' @examples 
#'     aves::extract_genus('Gallus gallus')
#'     aves::extract_genus('Gallus_gallus_2_TSA')
#' 
#' 
#' @export
extract_genus <- function(name) {
	stringr::str_extract(name, '^[:alpha:]+');
}


#' @title bind columns of taxonomy information
#'
#' @param data data.frame. must contain a column named `species`, whose first
#'   word is used to match genus name (in case of `'genus1 Sp'``). see
#'   [extract_genus]
#' @param taxonomy data.frame. taxonomy database to be used.
#'
#' @return data.frame
#'
#' @examples aves::bind_taxonomy(aves::genome)
#'
#' @export
bind_taxonomy <- function(data) {
	genus <- extract_genus(data$species);
	tax_subset <- aves::taxonomy %>% {.[match(genus, .$genus), 1:3]};
	dplyr::bind_cols(tax_subset, data) %>% dplyr::as_tibble();
}


#' @title correct species name error in NCBI genome & assembly database
#'
#' @param species character. species scientific name
#'
#' @return species
#' @export
#'
rectify_ncbi_species <- function(species) {
	species %<>% stringr::str_replace('Urile pelagicus', 'Phalacrocorax pelagicus');
    species %<>% stringr::str_replace('Nannopterum auritus', 'Phalacrocorax auritus');
    species %<>% stringr::str_replace('Nannopterum brasilianus', 'Phalacrocorax brasilianus');
	species %<>% stringr::str_replace('Nannopterum harrisi', 'Phalacrocorax harrisi');
	species %<>% stringr::str_replace('Emberiza melanocephala', 'Granativora melanocephala');
	species %<>% stringr::str_replace('Paradoxornis webbianus', 'Sinosuthora webbianus');
	species %<>% stringr::str_replace('Phylloscopus trochiloides', 'Seicercus trochiloides');
	species %<>% stringr::str_replace('Uraeginthus granatina', 'Granatina granatina');
	species %<>% stringr::str_replace('Gallirallus okinawae', 'Hypotaenidia okinawae');
	species %<>% stringr::str_replace('Phylloscopus plumbeitarsus', 'Seicercus plumbeitarsus');
	#species %<>% stringr::str_replace('', '');
	
	species;
}





#' @title show data distribution among orders
#'
#' @param data data.frame. must contain a column named `order`, such as `omics`.
#' 
#' @return `NULL`
#' 
#' @export
bin_by_order <- function(data) {
	data %>% {table(.$order)} %>% dplyr::as_tibble() %T>% 
		{print(dplyr::filter(., n > 2))} %>% 
		{print(dplyr::filter(., n > 3))};
	
	NULL;
}


#' @title preserve only the first row for same species
#'
#' @param data data.frame. must contain a column named `order`, such as `omics`.
#'
#' @return data.frame.
#' @export
#'
#' @examples aves::uniq_species(aves::rna)
# uniq_species <- function(data) {
# 	plyr::ddply(data, 'species', . %>% dplyr::slice(1));
# }


#' @title show orders significantly influenced (e.g., 2 data becomes 3 data) by a certain source
#'
#' @param data data.frame. must contain a column named `order`, such as `omics`.
#' @param Source string. such as `'SRA'` or `'genome'`
#' 
#' @return `NULL`
#' 
#' @export
diff_order <- function(data, Source) {
	order_3_no_SRA <- data %>% dplyr::filter(source != Source) %>% {table(.$order)} %>% 
	dplyr::as_tibble() %>% dplyr::filter(., n > 2) %>% {.[[1]]}

	order_3_SRA <- data %>% {table(.$order)} %>% 
		dplyr::as_tibble() %>% dplyr::filter(., n > 2) %>% {.[[1]]}
	
	
	order_4_no_SRA <- data %>% dplyr::filter(source != Source) %>% {table(.$order)} %>% 
		dplyr::as_tibble() %>% dplyr::filter(., n > 3) %>% {.[[1]]}
	
	order_4_SRA <- data %>% {table(.$order)} %>% 
		dplyr::as_tibble() %>% dplyr::filter(., n > 3) %>% {.[[1]]}
	
	return(union(setdiff(order_3_SRA, order_3_no_SRA), setdiff(order_4_SRA, order_4_no_SRA)))
}



#' @title summary data distribution using the YAML format
#'
#' @param data data.frame. must contain a column named `species`, such as
#'   `aves::protein`. would be transformed by [bind_taxonomy] if doesn't contain a
#'   column named `order`
#'
#' @return string. content of the yaml file. I suggest [readr::write_file()]
#' @export
#'
#' @examples aves::summary_as_yaml(aves::protein)
summary_as_yaml <- function(data) {
	bin <- function(df) {
		n_protein <- dplyr::filter(df, source %in% c('protein', 'genome')) %>% nrow()
		n_rna <- dplyr::filter(df, source %in% c('TSA', 'SRA')) %>% nrow()
		paste0(n_protein, 'protein + ', n_rna, 'rna')
	}

	data %<>% {if ('order' %in% colnames(.)) . else bind_taxonomy(.)};
	
	result <- data;
	result$order   %<>% {paste0(., ' -- ', plyr::daply(data, 'order', bin)[.])};
	result$family  %<>% {paste0(., ' -- ', plyr::daply(data, 'family', bin)[.])};
	result$genus   %<>% {paste0(., ' -- ', plyr::daply(data, 'genus', bin)[.])};
	result$species %<>% {paste0(., ' -- ', plyr::daply(data, 'species', bin)[.])};
	
	result %>% 
		plyr::dlply(
			'order', . %>% plyr::dlply(
				'family', . %>% plyr::dlply(
					'genus', . %>% plyr::dlply('species', . %>% {NULL})))) %>% 
		yaml::as.yaml(indent = 8)
}





