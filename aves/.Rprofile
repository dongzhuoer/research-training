if (file.exists('~/.Rprofile')) source('~/.Rprofile')

options(mcapomorphy.wd = '/path/to/data/aves') #" working directory on workstation
options(mcapomorphy.url = 'http://hostname:port') #" url for RStudio Server
options(mcapomorphy.max_iterate = 256L) #" url for RStudio Server




if (Sys.getenv('TRAVIS') == 'true' && Sys.getenv('CI') == 'true') {

} else {
    #doParallel::registerDoParallel(parallel::makeForkCluster());
    options(defaultPackages = c(options('defaultPackages')[[1]], 'tidyverse', 'xml2', 'rvest', 'parallel'));
}



