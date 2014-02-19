\name{getgeneinfo}
\alias{getgeneinfo}
\title{defines chromosome offsets for plotting multi chromosomal plot (eg plotManhattan)}
\usage{
getgeneinfo(biomart = "ensembl", host = NULL,
  dataset = "hsapiens_gene_ensembl", attributes = c("chromosome_name",
  "start_position", "end_position", "exon_chrom_start", "exon_chrom_end",
  "external_gene_id", "strand"), filters = c("chromosome_name", "start",
  "end", "with_ox_refseq_mrna"), values)
}
\arguments{
  \item{biomart}{biomart database to connect to}

  \item{host}{biomart host to connect to if different then
  www.biomart.org}

  \item{dataset}{biomart dataset you want to use}

  \item{attributes}{biomart attributes you want to
  retrieve.}

  \item{filters}{biomart filters (one or more) that should
  be used in the query.}

  \item{values}{values for the filters}
}
\description{
defines chromosome offsets for plotting multi chromosomal
plot (eg plotManhattan)
}

