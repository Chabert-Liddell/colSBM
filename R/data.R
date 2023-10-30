#' A collection of 8 food webs
#'
#' A dataset containing 8 stream food webs coming from the same study.
#'
#' Reference: Thompson, R. M. and Townsend, C. R. (2003),
#' Impacts on stream food webs of native and exotic forest:
#' an intercontinental comparison. Ecology 84(1): 145-161
#'
#' @format A named list of 8 binary matrices:
#' \describe{
#'   \item{M_Martins}{105 species and 343 feeding interactions}
#'   \item{NC_Cooper}{58 species and 126 feeding interactions}
#'   \item{NC_Herlzier}{71 species and 148 feeding interactions}
#'   \item{NZ_Venlaw}{69 species and 187 feeding interactions}
#'   \item{NZ_Berwick}{79 species and 240 feeding interactions}
#'   \item{NZ_North_Col}{78 species and 241 feeding interactions}
#'   \item{NZ_Powder}{78 species and 268 feeding interactions}
#'   \item{NZ_Trib}{98 species and 629 feeding interactions}
#' }
#' @source \url{http://www.web-of-life.es/}
"foodwebs"

#' A collection of 15 plant-pollinator bipartite networks
#'
#' A dataset of 15 plant-pollinator bipartite networks extracted from the same
#' study.
#'
#' Reference: Doré M., Fontaine C., Thébault E. (2021),
#' Relative effects of anthropogenic pressures, climate, and sampling design on
#' the structure of pollination networks at the global scale
#'
#' @format A named list of 15 binary adjacency matrices:
#' \describe{
#'  \item{arroyo1982_1+arroyo1982_2+arroyo3}{113 pollinator species, 131 plant species and 544 pollination interactions}
#'  \item{dupont2003}{38 pollinator species, 11 plant species and 106 pollination interactions}
#'  \item{eberling1999}{118 pollinator species, 23 plant species and 238 pollination interactions}
#'  \item{herrera1988}{179 pollinator species, 26 plant species and 412 pollination interactions}
#'  \item{inouye1988}{88 pollinator species, 41 plant species and 279 pollination interactions}
#'  \item{kato1990}{679 pollinator species, 93 plant species and 1207 pollination interactions}
#'  \item{medan2002ld}{45 pollinator species, 21 plant species and 83 pollination interactions}
#'  \item{medan2002rb}{72 pollinator species, 23 plant species and 125 pollination interactions}
#'  \item{olensen2002aig}{13 pollinator species, 14 plant species and 52 pollination interactions}
#'  \item{olensen2002flo}{12 pollinator species, 10 plant species and 30 pollination interactions}
#'  \item{ramirez1992}{53 pollinator species, 28 plant species and 109 pollination interactions}
#'  \item{small1976}{34 pollinator species, 13 plant species and 141 pollination interactions}
#'  \item{vazquez2002}{90 pollinator species, 14 plant species and 164 pollination interactions}
#'  \item{petanidou1991}{666 pollinator species, 131 plant species and 2933 pollination interactions}
#'  \item{ramirez1989}{49 pollinator species, 48 plant species and 156 pollination interactions}
#' }
#' @source \url{https://zenodo.org/records/4300427}
"dorebipartite"
