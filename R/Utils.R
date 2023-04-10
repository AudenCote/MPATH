require(readxl)


read_sheet <- function(path, xls_sheet = 'NA'){
  ext <- file_ext(path)
  if(ext == 'tsv'){
    data.frame(read_tsv(path))
  }else if(ext == 'csv'){
    data.frame(read_csv(path))
  }else if(ext == 'xls' | ext == 'xlsx'){
    if(xls_sheet != 'NA'){
      data.frame(read_excel(path, sheet = xls_sheet))
    } else{
      data.frame(read_excel(path))
    }
  }
  else{
    logger('File format not recognized. Please input a tsv, csv, or .xlsx file', 'error')
  }
}
