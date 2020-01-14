## parse hamstr output and generate fasta per gene


1. hamstr-out
1. species (use `in_tsa` of `sra` and `\d{9}` of dna to check duplication)
1. mine


```r
#options(mc.cores = 256L);
#c('aln', 'mafft', 'muscle', 'aliscore') %>% paste0('gene/', .) %>% 
c('hamstr-out', 'species', 'gene') %>% paste0('paup/', .) %>% 
    plyr::l_ply(. %>% {if(dir.exists(.)) unlink(., recursive = T); dir.create(., recursive = T)})

# merge hamstr output
str_extract(dir('hamstr', '\\.out'), '(?<=_)\\w+(?=_aves)') %T>% 
    mclapply(. %>% {
        input <- c(
            paste0('hamstr', '/hamstrsearch_', ., '_aves.out'),
            paste0('hamstr-ext', '/hamstrsearch_', ., '_', ., '.out')
        );
        output <- paste0('apocode/hamstr-out/', ., '.faa');

        mcapomorphy::hamstr_out_to_fasta(input, output);
    }, mc.cores = 256L) -> dev.null;

# select best data per species for mine
dir('apocode/hamstr-out', full.names = T) %>% 
    mclapply(. %>% {
        data_frame(src = ., n = read_lines(.) %>% {length(.)/2});
    }) %>%
    bind_rows() %>% 
    mutate(species = str_extract(src, '[:alpha:]+_[:alpha:]+')) %>%
	plyr::dlply('species', . %>% {filter(., n == max(.$n))}) %>% 
    mclapply(. %>% {
        dest <- paste0('apocode/species/', .$species, '.faa');
        file.copy(.$src, dest, overwrite = T);
    }, mc.cores = 256L) -> dev.null;

# organize sequences into OG
dir('apocode/species', full.names = T) %>% 
    {.[!str_detect(., 'Pygoscelis_papua')]} %>%
    {.[!str_detect(., 'Otus_bakkamoena')]} %>%
    mclapply(. %>% {
        species = basename(.) %>% str_extract('\\w+');
        fasta <- bioinfor::read_fasta(.);
        data_frame(species, gene = names(fasta), sequence = fasta);
    }, mc.cores = 256L) %>% 
    bind_rows() %>%
    plyr::dlply('gene') %>%
    mclapply(. %>% {
        fasta <- str_replace_all(.$sequence, '\\*', 'X') ; #" muscle and paup both don't recognize '*' 
        names(fasta) <- .$species;
        
        path <- paste0('apocode/gene/', .$gene[1], '.aln');
        bioinfor::write_fasta(fasta, path);
    }, mc.cores = 256L) -> dev.null;

# check whether sequence header in aln is duplicated
dir('apocode/gene', pattern = 'aln$', full.names = T) %>%
    mclapply(biozhuoer::read_fasta) %>%
    plyr::llply(. %>% {.$name} %>% {.[duplicated(.)]}) %>% plyr::laply(length)  %>% {. != 0} %>% which
# .Last.value -> bad
# dir('mine/gene', pattern = 'aln$', full.names = T)[bad]
```



## for verify

```r
dir('paup/hamstr-out', full.names = T) %>% 
    mclapply(. %>% {
        species = basename(.) %>% str_extract('\\w+');
        fasta <- bioinfor::read_fasta(.);
        data_frame(species, gene = names(fasta), sequence = fasta);
    }, mc.cores = 256L) %>% 
    bind_rows() %>%
    plyr::dlply('gene') %>%
    mclapply(. %>% {
        fasta <- str_replace_all(.$sequence, '\\*', 'X') ; #" muscle and paup both don't recognize '*' 
        names(fasta) <- .$species;
        
        path <- paste0('verify/gene/', .$gene[1], '.aln');
        bioinfor::write_fasta(fasta, path);
    }, mc.cores = 256L) -> dev.null;
```









## align and trim

```r
mclapply(
    dir('apocode/gene', 'aln$', full.names = T), 
    . %>% {
        dest <- str_replace(., 'aln$', 'mafft');
        system2('mafft', c('--auto', '--amino', '--anysymbol', .), dest, F);
        if(file.size(dest) == 0L) bioinfor::read_fasta(.) %>% bioinfor::write_fasta(dest);    
    }
) -> dev.null;

mclapply(
    dir('apocode/gene', 'mafft$', full.names = T), 
    . %>% {
        dest <- str_replace(., 'mafft$', 'muscle');
        system2('muscle',c('-refine', '-in', .), dest, F);       
        #if(file.size(dest) == 0L) bioinfor::read_fasta(.) %>% bioinfor::write_fasta(dest);    
    }#, mc.cores = as.integer(getOption('mc.cores')/0.8)
) -> dev.null;

# mclapply(
#     dir('apocode/gene', 'muscle$', full.names = T), 
#     . %>% {bioinfor::aliscore(c('-i', ., '-r 6670'), F)} #" C(116,2) = 116*115/2 = 6670
# ) -> dev.null;
# #" I remove empty List files so that every List file shall have a corresponding ALICUT_ file
# dir('apocode/gene', 'List|Profile', full.names = T) %>% {file.remove(.[file.size(.) <= 2L])}

# biozhuoer::alicut('apocode/gene', '-s');

# mclapply(
#     dir('apocode/gene', '^E\\w+.muscle$', full.names = T), 
#     . %>% {
#         dest <- str_replace(., 'muscle', 'trim');
#         alicut_out <- str_replace(., 'gene/', 'gene/ALICUT_');

#         if(file.exists(alicut_out)) {
#             file.copy(., dest)
#         } else {
#             bioinfor::read_fasta(.) %>% bioinfor::write_fasta(dest);
#         }
#     } 
# ) -> dev.null;


dir('apocode/gene',full.names = T) %>% {data_frame(name = ., size = file.size(.))} %>% slice(., order(.$size)) %T>% print(n=100)
dir('apocode/gene', 'aln$') %>% length == dir('apocode/gene', 'muscle$') %>% length
# dir('apocode/gene', 'List') %>% length == dir('apocode/gene', 'ALICUT_') %>% length - 1
# dir('apocode/gene', 'aln$') %>% length == dir('apocode/gene', 'trim$') %>% length
```


## concatenate

```r
#we only preserve OG which is contained by above 90% species
all_species <- dir('paup/species') %>% str_replace('.faa', '') %>% sort;
dir('apocode/gene', pattern = 'trim$', full.names = T) %>%
    mclapply(bioinfor::read_fasta) %>% 
    {.[mclapply(., length) > dir('paup/species') %>% length * 0.75]} %>% 
    mclapply(. %>% {
        # construct missing species using 'gapped' sequence
        missed_species <- names(.) %>% setdiff(all_species, .);
        .[missed_species] = .[1] %>% nchar %>% str_dup('-', .) #" fasta is already aligned, the first one is enough to represent all sequence
        . = .[order(names(.))]
    }, mc.cores = 256) %>% 
    do.call(paste0, .) %>%    
    {names(.) = all_species;.} %>% 
    bioinfor::write_fasta('paup/aves.fasta');





dir('apocode/gene', pattern = 'trim$', full.names = T) %>%
    mclapply(biozhuoer::read_fasta) %>% 
    mclapply(y[4305], . %>% {
        # supplement sequence for missing species using 'gapped' sequence
        all_species <- dir('paup/species') %>% str_replace('.faa', '') %>% sort;
        #" fasta is already aligned, the first one is enough to represent all sequence
        forge_seq <- .$seq[1] %>% nchar %>% str_dup('-', .);
        missed <- tibble(name = setdiff(all_species, .$name), seq = forge_seq)
        result <- if (nrow(.) < length(all_species)) bind_rows(., missed) else .
        result %<>% slice(., order(.$name))
    })->z %>%
    {tibble(name = .[[1]]$name, seq = do.call(paste0, mclapply(., . %>% {.$seq})))} %>% 
    biozhuoer::write_fasta('paup/aves.matrix')
```



## paup
```r
fasta_sub <- function(fasta, ...) {
    print(names(fasta))
    result <- str_sub(fasta, ...);
    names(result) <- names(fasta);
    result
}
bioinfor::read_fasta('paup/aves.fasta') %>% fasta_sub(1001, 2000) %>% bioinfor::write_fasta('paup/aves1000.fasta')


mcapomorphy::write_paup('paup/aves.fasta', aves::pkg_file('omics.tre'), 
                 aves::apomorphy_mine('aves.paup', aves::outgroup$species), 
                 'paup/aves.nexus')

system.time(bioinfor::paup(c('paup/aves.nexus', '-n', '-u'), stdout = F));



#system2('iconv', '-f US-ASCII -t UTF-8 paup/aves.paup > paup/aves.paup.utf8');
```

`paup4b10 paup/aves.nexus -n -u`



