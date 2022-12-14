% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_RNASeq.R
\name{run_DESeq2}
\alias{run_DESeq2}
\title{Title A pipeline of run RNA-seq DEGs analysis and mini database of Schistosoma japonicum}
\usage{
run_DESeq2(
  count_data,
  groupList = c("GFP", "dsRNA"),
  group_a_name = NULL,
  group_a_samples = NULL,
  group_b_name = NULL,
  group_b_samples = NULL,
  group_c_name = NULL,
  group_c_samples = NULL,
  group_d_name = NULL,
  group_d_samples = NULL,
  group_e_name = NULL,
  group_e_samples = NULL,
  group_f_name = NULL,
  group_f_samples = NULL,
  group_g_name = NULL,
  group_g_samples = NULL,
  group_h_name = NULL,
  group_h_samples = NULL,
  pheno_file = "",
  batch = NULL
)
}
\arguments{
\item{count_data}{file path of count matrix which generated by previous RNA-seq pipeline,
eg: "2.Results/SnakePipe_outputs_test/featureCounts/counts.tsv"}

\item{groupList}{}

\item{group_a_name}{group a name eg:GFP}

\item{group_a_samples}{a list of samples in group A you choose, must in the colname of count_data}

\item{group_b_name}{group b name eg:dsNFY}

\item{group_b_samples}{a list of samples in group B you choose, must in the colname of count_data}

\item{group_c_name}{group names}

\item{group_c_samples}{a list of samples in group , must in the colname of count_data}

\item{group_d_name}{group names}

\item{group_d_samples}{a list of samples in group , must in the colname of count_data}

\item{group_e_name}{group names}

\item{group_e_samples}{a list of samples in group , must in the colname of count_data}

\item{group_f_name}{group names}

\item{group_f_samples}{a list of samples in group , must in the colname of count_data}

\item{group_g_name}{group names}

\item{group_g_samples}{a list of samples in group , must in the colname of count_data}

\item{group_h_name}{group names}

\item{group_h_samples}{a list of samples in group , must in the colname of count_data}

\item{pheno_file}{if group > 8, please provide pheno info in a csv file format.}

\item{batch}{if have batch effect, eg c(1,1,1,2,2,1),1~ batch1,2~ batch2}
}
\value{
a dds class for downstream analysis
}
\description{
Title A pipeline of run RNA-seq DEGs analysis and mini database of Schistosoma japonicum
}
\examples{
dds <- run_DESeq2(
  count_data = "~/Projects/RNASeq_CGW_20221014/2.Results/SnakePipe_outputs/featureCounts/counts.tsv",
  group_a_name = 'GFP',
  group_a_samples  = c("G.F1","G.F3","G.M"),
  group_b_name = 'WSS',
  group_b_samples = c("W.F1","W.F2","W.M1")
)
}
