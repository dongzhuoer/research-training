
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE);
```

## prerequisite

```
hamstr-data
  |---- aves
  |---- blast_dir
  |---- ext
hamstr
hamstr-ext
```


## data


```r
# protein and outgroup
mclapply(
    dir('protein'),
    . %>% {
        input  <- paste0('protein/', .);
        output <- paste0('hamstr/', .) %>% str_replace('.gz$', '');
        mcapomorphy::hamstr_from_ncbi_protein(input, output);
    }
) -> dev.null



# genome
mclapply(
    dir('augustus'),
    . %>% {
        input  <- paste0('augustus/', .);
        output <- paste0('hamstr/', .) %>% str_replace(., '.gff$', '.faa');
        mcapomorphy::hamstr_from_augustus(input, output);
    }
) -> dev.null



# TSA

# go to aves/process.Rmd to prepare input
# .gz or .fna

mclapply(
    dir('transcript'),
    . %>% {
        input  <- paste0('transcript/', .);
        output <- paste0('hamstr/', .) #%>% str_replace('.gz$', '');
        mcapomorphy::hamstr_from_TSA(input, output);
    }
) -> dev.null
       


# SRA
filter(aves::sra, id %in% str_extract(dir('oases/output'), '\\w+'))

mclapply(
    dir('augustus'),
    . %>% {
        input  <- paste0('oases/output/', .);
        output <- paste0('hamstr/', .) %>% str_replace(., '.fa$', '.faa');
        mcapomorphy::hamstr_from_augustus(input, output);
    }
) -> dev.null
```


## run hamstr

```r
str_extract(dir('hamstr', '\\.f[an]a$'), '\\w+') %>% 
    {setdiff(., str_extract(dir('hamstr', '\\.out'), '(?<=_)\\w+(?=_aves)'))} %>%
    aves::hamstr_command() %>% 
    c('#!/bin/bash', .) %>% write_lines('hamstr.sh')
```

```bash
cp hamstr/*aves.out hamstr/output_top10
```

> run time profile

```yaml
protein: 2~39
gen{ome,ine}: 11~46
outgroup: 15~118(human)
TSA: 28~36,258
```



## extension (other than top 10 hmmer hit)

```r
system('rm -r hamstr/hmm_dir_*');
system('rm -r hamstr/fa_dir_*');

mclapply(
    dir('hamstr', 'aves.out', full.names = T) ,
    . %>% {
        name <- str_extract(., '(?<=_)\\w+(?=_aves)');
        source_dir <- 'hamstr-data/aves/hmm_dir/';
        dest_dir <- paste0('hamstr/hmm_dir_', name, '/aves/hmm_dir/') %T>% dir.create(recursive = T);
        
        todo <- read_lines(.) %>% str_extract('^\\w+') %>% setdiff(dir(source_dir), .) %>% c('../aves.fa')
        
        file.copy(paste0(source_dir, todo), paste0(dest_dir, todo));
    }, mc.cores = 128
) -> dev.null

dir('hamstr', '^hmm_dir', full.names = T) %>% paste0('/aves/hmm_dir') %>% plyr::laply(. %>% dir %>% length) %>% sort
```

```bash
cp hamstr/*aves.out hamstr/output_all
```

## check whether fasta file contains duplicated header

```{r, eval=FALSE}
dir('protein10/', 'faa$', full.names = T) %>% plyr::laply(.parallel = T,. %>% read_lines %>% str_subset('^>') %>% duplicated %>% which %>% length)
dir('TSA10/', 'fna$', full.names = T) %>% plyr::laply(.parallel = T,. %>% read_lines %>% str_subset('^>') %>% duplicated %>% which %>% length)
```

dir('genome10', 'out', full.names = T) %>% {na.omit(str_extract(., '(?<=search_)\\w+(?=_aves)'))} %in% 
str_replace(filter(aves::pool, source == 'genein')$name, 'protein', 'genome')}






```r
aves::sra$id %>% str_replace_all('[_\\.]', '')

```