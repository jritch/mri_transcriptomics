writeTableWrapper <- function (my_table,dir,filename,col_names=NULL,rows=NULL,col_indices=NULL) {
  if (length(rows) > 0){
    my_table = head(my_table,n=rows)
  }
  if (length(col_indices) > 0) {
    my_table=my_table[col_indices]
  }
  if (length(col_names) > 0) {
   names(my_table) <- col_names
   write.csv(my_table, file.path(dir,filename))
  }

  write.csv(my_table, file.path(dir,filename))
  return(my_table)
}