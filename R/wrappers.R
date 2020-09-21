#' @export
runOcularTxomeVis<- function() {
    appDir <- system.file("app",  package = "OcularTxome")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
    } else if(!file.exists(paste0(appdir, '/app_data/shiny_data.Rdata'))){
        stop("Please Download app data. Run OcularTxome::downloadAppData() ")
    }

    shiny::runApp(appDir, display.mode = "normal")
}
#' @export
downloadAppData <- function(){
    appDir <- system.file("app",  package = "OcularTxome")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
    }
    download.file(url ='http://hpc.nih.gov/~mcgaugheyd/ocular_transcriptomes_shiny/app_data.tar.gz',
                destfile = paste0(appDir, '/app_data.tar.gz'),
                quiet = T)
    untar(paste0(appDir, '/app_data.tar.gz'), exdir = appDir)
    unlink(paste0(appDir, '/app_data.tar.gz'))
}
#' @export
downloadAnnotation <- function(path=''){
    download.file(url = 'http://hpc.nih.gov/~mcgaugheyd/ocular_transcriptomes_shiny/dl_data.tar.gz',
                  dest = paste0(path, 'dl_data/tar.gz') )
    untar( paste0(path, 'dl_data/tar.gz') )
    unlink( paste0(path, 'dl_data/tar.gz') )

}
