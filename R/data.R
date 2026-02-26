#' Spatial Transcriptomics Data for the Medial Prefrontal Cortex (mPFC) of Three Mice
#'
#' A publicly available dataset from the spatial transcriptomics technology STARMAP, which consists of three tissue sections obtained from the medial prefrontal cortex of the brain from three different mice.
#'
#' @format A list of two objects:
#' \describe{
#' \item{starmap_cnts}{a list of three gene expression count matrices (dimensions: number of genes X number of cells), one for each tissue sample}
#' \item{starmap_info}{a list of metadata for each tissue sample (dimensions: number of cells X 4), containing, in order, the x- and y-coordinates of the cell centroids, and the manually annotated cell type and spatial domain label for each cell}
#' }
#' @references
#' X. Wang, W. E. Allen, M. A. Wright, E. L. Sylwestrak, N. Samusik, S. Vesuna, K. Evans, C. Liu, C. Ramakrishnan, J. Liu, G. P. Nolan, F.-A. Bava, and K. Deisseroth. Three-dimensional intact-tissue sequencing of single-cell transcriptional states. Science, 361(6400):eaat5691, July 2018. doi: 10.1126/science.aat5691.
#' Z. Li and X. Zhou. BASS: multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies. Genome Biology, 23(1):168, Aug. 2022. ISSN 1474-760X. doi: 10.1186/s13059-022-02734-7.
"starmap_mpfc"



# TODO: make dataset consisting of sliver of 2 tissue samples from original data
