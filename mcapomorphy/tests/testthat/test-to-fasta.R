testthat::context("Testing to-fasta.R")
setwd(here::here(''))  # workspace is reset per file




# test_that("Testing augustus_gff_to_fasta", {
# 	augustus_demo <- tempfile() %T>% augustus_gff_to_fasta('data-raw/augustus-demo.gff', .) %>% biozhuoer::read_fasta();

# 	expect_true(identical(augustus_demo$name, c('1', '2')));
# 	expect_true(identical(nchar(augustus_demo$seq), c(604L, 167L)));
# });



testthat::test_that("Testing hamstr_out_to_fasta", {
    # I make line 3 duplicate of line 2 to test `hamstr_out_to_fasta`
	chicken <- tempfile() %T>% hamstr_out_to_fasta('inst/extdata/chicken.out', .) %>% biozhuoer::read_fasta();

	testthat::expect_true(identical(chicken$name, c('EOG090F0013', 'EOG090F0028')));
	testthat::expect_true(identical(nchar(chicken$seq), c(3478L, 3397L)));
});
