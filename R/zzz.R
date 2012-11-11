setClass("prediction",
         representation(predictions = "list",
                        labels      = "list",
                        cutoffs     = "list",
                        fp          = "list",
                        tp          = "list",
                        tn          = "list",
                        fn          = "list",
                        n.pos       = "list",
                        n.neg       = "list",
                        n.pos.pred  = "list",
                        n.neg.pred  = "list"))

setClass("performance",
         representation(x.name       = "character",
                        y.name       = "character",
                        alpha.name   = "character",
                        x.values     = "list",
                        y.values     = "list",
                        alpha.values = "list" ))

#setMethod("plot",signature(x="performance",y="missing"),
#          function(x,y,...) {
#              .plot.performance(x,...)
#          })
