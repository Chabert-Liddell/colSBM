library(profvis)
devtools::load_all()
data(dorebipartite)

netlist <- dorebipartite[1:4]
p <- profvis({
    fit <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "iid",
        net_id = names(netlist),
        nb_run = 3L,
        global_opts = list(backend = "no_mc")
    )
})
filepath <- file.path(paste0(tempfile(), ".html"))
htmlwidgets::saveWidget(p, filepath)
browseURL(filepath)
