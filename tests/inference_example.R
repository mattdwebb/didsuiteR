if (requireNamespace("fasttwfeR", quietly = TRUE)) {
  source(system.file("examples", "inference.R", package = "didsuiteR"))
} else {
  message("Inference example requires the optional fasttwfeR package; release integration installs it.")
}
