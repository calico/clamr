#' Get SQL database schema
#'
#' @param sql_con a DBI connection to any sql database
#'
#' @export
sql_get_schema <- function(sql_con) {
  db_tables <- DBI::dbListTables(sql_con)

  output_list <- list()
  for (a_db_table in db_tables) {
    db_nrow <- dplyr::tbl(
      sql_con,
      dbplyr::sql(paste0("SELECT Count(*) FROM ", a_db_table))
    ) %>%
      dplyr::collect() %>%
      unlist() %>%
      unname()

    if (db_nrow == 0) {
      # no content in this table
      next
    }

    db_head <- suppressWarnings(
      dplyr::tbl(
        sql_con,
        dbplyr::sql(paste0("SELECT * FROM ", a_db_table))
      ) %>%
        dplyr::collect(n = 10)
    )

    output_list[[a_db_table]] <- list()
    output_list[[a_db_table]]$db_nrow <- db_nrow
    output_list[[a_db_table]]$db_head <- db_head
  }

  output_list <- output_list %>%
    purrr::transpose()

  output_list$db_nrow <- output_list$db_nrow %>%
    {
      tibble::tibble(table = names(.), nrows = unlist(.))
    }

  output_list
}

#' Safe SQLite Con
#'
#' Connect to an SQLite database.
#'
#' @param sqlite_path path to an .sqlite database
#' @param unlock if TRUE then unlock database if database is currently locked
#'
#' @return \code{sqlite_con}: a connection to the mzroll/sqlite database
#'
#' @export
create_sqlite_con <- function(sqlite_path, unlock = TRUE) {
  checkmate::assertFileExists(sqlite_path)
  checkmate::assertLogical(unlock, len = 1)

  sqlite_con <- try(
    DBI::dbConnect(RSQLite::SQLite(), sqlite_path, synchronous = NULL),
    silent = TRUE
  )

  if ("try-error" %in% class(sqlite_con)) {
    warning(
      "SQLite database cannot be used: database is unreadable, ",
      attributes(sqlite_con)$condition$message,
      "\n"
    )
    return(NULL)
  }

  # test whether sql is locked and unlock if needed
  is_sqlite_lock <- test_sqlite_lock(sqlite_con)

  if (is_sqlite_lock) {
    if (unlock) {
      sqlite_con <- unlock_sqlite(sqlite_con)
    } else {
      stop(
        "Database is locked,
        call create_sqlite_con() with unlock = TRUE to unlock"
      )
    }
  }

  return(sqlite_con)
}

test_sqlite_lock <- function(sqlite_con) {
  if (!("SQLiteConnection" %in% class(sqlite_con))) {
    stop("\"sqlite_con\" must be an sqlite database connection")
  }

  try_unlock <- try(DBI::dbListTables(sqlite_con), silent = TRUE)

  if (
    "try-error" %in% class(try_unlock) &&
      stringr::str_detect(try_unlock[1], "database is locked")
  ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Unlock a locked sqlite database
#'
#' This is a total hack, where I copy the database to /tmp and then over-write
#'   the original file
#'
#' @param sqlite_con a DBI connection to any a sqlite database
#'
#' @return an unlocked sqlite_con
unlock_sqlite <- function(sqlite_con) {
  if (!("SQLiteConnection" %in% class(sqlite_con))) {
    stop("\"sqlite_con\" must be an sqlite database connection")
  }

  locked_db_path <- sqlite_con@dbname
  if (file.access(locked_db_path, mode = 2) == -1) {
    stop("Database is not writeable and cannot be automatically unlocked")
  }

  tmp_basename <- glue::glue("tmb_db_{floor(runif(1, 0, 1000000))}.sqlite")
  tmp_file <- file.path(dirname(locked_db_path), tmp_basename)
  if (file.access(tmp_file, mode = 2) == -1) {
    # try tmp directory (but may run out of memory)
    tmp_file <- file.path("/tmp", tmp_basename)
  }

  warning(
    "Unlocking locked database and recreating connection: ",
    locked_db_path
  )

  DBI::dbDisconnect(sqlite_con)
  system(glue::glue("cp {locked_db_path} {tmp_file}"))
  system(glue::glue("rm {locked_db_path}"))
  system(glue::glue("mv {tmp_file} {locked_db_path}"))
  sqlite_con <- DBI::dbConnect(
    RSQLite::SQLite(),
    locked_db_path,
    synchronous = NULL
  )

  return(sqlite_con)
}

#' SQLite Write
#'
#' Write to an sqlite database testing for and unlocking database if needed
#'
#' @inheritParams unlock_sqlite
#' @param name name of table to write
#' @param value values of table
#' @param ... additionally arguments to pass to \link[DBI]{dbWriteTable}.
#'
#' @return returns 0 invisibly and writes to database as a side effect
#'
#' @export
sqlite_write <- function(sqlite_con, name, value, ...) {
  if (test_sqlite_lock(sqlite_con)) {
    sqlite_con <- create_sqlite_con(sqlite_con@dbname)
  }

  DBI::dbWriteTable(sqlite_con, name, value, ...)

  if (test_sqlite_lock(sqlite_con)) {
    sqlite_con <- create_sqlite_con(sqlite_con@dbname)
  }

  return(invisible(0))
}

#' SQLite Copy
#'
#' Copy an SQLite database from one location to another, testing for
#'   malformity.
#'
#' @param sqlite_from an sqlite-con created by \link{create_sqlite_con}
#' @param sqlite_to an sqlite-con created by \link{create_sqlite_con}
#' @param overwrite overwrite an existing \code{sqlite_to}
#'
#' @return 0 invisibly if copying occurred, 1 otherwise
#'
#' @export
sqlite_copy <- function(sqlite_from, sqlite_to, overwrite = FALSE) {
  checkmate::assertFileExists(sqlite_from)
  checkmate::assertLogical(overwrite, len = 1)

  invalid_source <- sqlite_malformed(sqlite_from)
  if (invalid_source) {
    stop(
      sqlite_from,
      " is malformed, please correct before calling sqlite_copy()"
    )
  }

  if (!overwrite) {
    if (file.exists(sqlite_to)) {
      invalid_dest <- sqlite_malformed(sqlite_to)

      if (!invalid_dest) {
        return(invisible(1))
      }
    }
  }

  file.copy(sqlite_from, sqlite_to, overwrite = TRUE)

  invalid_output <- sqlite_malformed(sqlite_to)
  if (invalid_output) {
    stop(sqlite_to, " output was malformed, try again")
  }

  return(invisible(0))
}

sqlite_malformed <- function(sqlite_path) {
  sqlite_con <- try(
    create_sqlite_con(sqlite_path, unlock = FALSE),
    silent = TRUE
  )

  if ("try-error" %in% class(sqlite_con)) {
    error_message <- attr(sqlite_con, "condition")[1]$message

    # locking is okay
    if (stringr::str_detect(error_message, "^Database is locked")) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  is_malformed <- try(DBI::dbListTables(sqlite_con), silent = TRUE) %>%
    {
      "try-error" %in% class(.)
    }
  DBI::dbDisconnect(sqlite_con)

  return(is_malformed)
}
