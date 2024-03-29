% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqlite.R
\name{create_sqlite_con}
\alias{create_sqlite_con}
\title{Safe SQLite Con}
\usage{
create_sqlite_con(sqlite_path, unlock = TRUE)
}
\arguments{
\item{sqlite_path}{path to an .sqlite database}

\item{unlock}{if TRUE then unlock database if database is currently locked}
}
\value{
\code{sqlite_con}: a connection to the mzroll/sqlite database
}
\description{
Connect to an SQLite database.
}
