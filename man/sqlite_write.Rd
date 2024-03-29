% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqlite.R
\name{sqlite_write}
\alias{sqlite_write}
\title{SQLite Write}
\usage{
sqlite_write(sqlite_con, name, value, ...)
}
\arguments{
\item{sqlite_con}{a DBI connection to any a sqlite database}

\item{name}{name of table to write}

\item{value}{values of table}

\item{...}{additionally arguments to pass to \link[DBI]{dbWriteTable}.}
}
\value{
returns 0 invisibly and writes to database as a side effect
}
\description{
Write to an sqlite database testing for and unlocking database if needed
}
