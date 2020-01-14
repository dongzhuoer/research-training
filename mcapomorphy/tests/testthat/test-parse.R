testthat::context("Testing parse.R")
setwd(here::here(''))  # workspace is reset per file


testthat::test_that("Testing parse_genome", {
    parse_genome_55342 <- tibble::data_frame(
        species = "Nannopterum harrisi", number = 1L,
        dna = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/173/475/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna.gz",
        rna = "", protein = ""
    );
    testthat::expect_true(identical(parse_genome('inst/extdata/genome55342.html'), parse_genome_55342));

    parse_genome_2793 <- tibble::data_frame(
        species = "Anas platyrhynchos", number = 2L,
        dna = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/355/885/GCF_000355885.1_BGI_duck_1.0/GCF_000355885.1_BGI_duck_1.0_genomic.fna.gz",
        rna = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/355/885/GCF_000355885.1_BGI_duck_1.0/GCF_000355885.1_BGI_duck_1.0_rna.fna.gz",
        protein = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/355/885/GCF_000355885.1_BGI_duck_1.0/GCF_000355885.1_BGI_duck_1.0_protein.faa.gz"
    );
    testthat::expect_true(identical(parse_genome('inst/extdata/genome2793.html'), parse_genome_2793));
});



test_that("Testing parse_avibase", {
    parse_avibase_demo <- tibble::tibble(
        order = "Coliiformes", family = "Coliidae", genus = "Colius",
        species = "Colius striatus", id = 'Avibase-1FDDABDB0D4421F9.html'
    );
    expect_true(identical(parse_avibase('inst/extdata/Avibase-1FDDABDB0D4421F9.html'), parse_avibase_demo));
});
