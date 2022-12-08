require(readxl)

logger <- function(text, type){

  if(type == 'error'){

  }else if(type == 'warning'){

  }else if(type == 'message'){

  }

}


read_sheet <- function(path, xls_sheet = 'NA'){
  ext <- file_ext(path)
  if(ext == 'tsv'){
    data.frame(read_tsv(path))
  }else if(ext == 'csv'){
    data.frame(read_csv(path))
  }else if(ext == 'xls' | ext == 'xlsx'){
    data.frame(read_excel(path, sheet = xls_sheet))
  }
  else{
    logger('File format not recognized. Please input a tsv, csv, or .xlsx file', 'error')
  }
}
