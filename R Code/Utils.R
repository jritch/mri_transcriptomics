writeTableWrapper <- function (table,dir,filename,col_names=NULL,rows=NULL,col_indices=NULL) {
  if (rows){
    table = head(table,n=rows)
  }
  if (col_indices) {
    table=table[col_indices]
  }
  if (length(col_names) > 0) {
   names(table) <- col_names
   write.csv(table, file.path(dir,filename))
  }

  write.csv(table, file.path(dir,filename))
}