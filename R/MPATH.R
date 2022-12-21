require(tidyverse)
require(tools)
require(cluster)
require(ggfortify)
require(factoextra)
require(NbClust)
require(purrr)
require(ggrepel)
require(jsonlite)
require(httr)
require(RColorBrewer)
require(grid)
require(gridExtra)
require(heatmaply)

source('Utils.R')

MPATH_Pipeline <- setRefClass('MPATH_Pipeline',

  fields = list(
    expression_file = 'character',
    expression_data = 'data.frame',
    mitocarta_data = 'list',
    log2fc_pval_dataframe = 'data.frame',
    volcano_plots = 'list',
    regulation_barplot = 'list',
    goi_expression = 'matrix',
    pca = 'list',
    pathways = 'list',
    test_hmdb = 'matrix',
    top_pathway_genes_l2fcp = 'data.frame'
  ),

  methods = list(

  initialize = function(...,
                          mitocarta_data = list(
                           'path' = '../Data/Human.MitoCarta3.0.xls',
                           'gene.sheet' = 'A Human MitoCarta3.0',
                           'gene.col' = 'EnsemblGeneID_mapping_version_20200130',
                           'pathway.sheet' = 'C MitoPathways'
                          ),
                          pathways = list(
                            'mitocarta_frequencies' = data.frame(),
                            'mitocarta_fisher_tests' = data.frame(),
                            'panther_frequencies' = data.frame()
                          )
                         ){
    callSuper(..., mitocarta_data = mitocarta_data, pathways = pathways)

    expression_data <<- data.frame(read_tsv(.self$expression_file)) %>%
     sample_n(1000)

    .self$expression_data[.self$expression_data == 0] <- NA

    expression_data <<- .self$expression_data %>%
      na.omit() %>%
      pivot_longer(!Gene, names_to = 'Sample', values_to = 'Exp') %>%
      separate(Sample, c('Sample', 'Replicate'), sep = '_') %>%
      filter(Gene %in% read_sheet(.self$mitocarta_data$path, .self$mitocarta_data$gene.sheet)[[.self$mitocarta_data$gene.col]]) %>%
      merge(.self$ensembl_query(.self$expression_data$Gene), by = 'Gene')
  },

  ensembl_query = function(ids){
    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    ids_json <- toJSON(ids)

    body <- list(api=1, ids=ids_json)
    r <- POST(url, body = body)

    output <- fromJSON( content(r, "text"), flatten=TRUE)
    output <- output[-which(sapply(output, is.null))]
    #fromJSON( content(r, "text"), flatten=TRUE)
    data.frame(Gene = names(output), Symbol = unname(unlist(output)))
  },


  panther_query = function(goi_list){
    url = "http://www.pantherdb.org/services/oai/pantherdb/enrich/overrep?"

    body <- list(
      'geneInputList' = paste(goi_list, collapse=","),
      'organism' = '9606',
      'annotDataSet' = 'ANNOT_TYPE_ID_PANTHER_PATHWAY',
      'enrichmentTestType' = 'FISHER',
      'correction' = 'FDR'
    )

    r <- POST(url, body = body, encode = "form", verbose())
    data.frame(fromJSON( content(r, "text"), flatten=TRUE)$results$result) %>%
      select(fold_enrichment | term.label) %>%
      rename(Pathway = term.label) %>%
      rename(Fold.Enrichment = fold_enrichment)
  },

  Volcano = function(benchmark_sample = NULL){
    if(is.null(benchmark_sample)){
      benchmark_sample <- unique(.self$expression_data$Sample)[1]
    }

    d1 <- filter(.self$expression_data, Sample == benchmark_sample) %>% mutate(Exp1 = Exp) %>% select(!Sample & !Exp)
    compared_samples <- unique(.self$expression_data$Sample)[unique(.self$expression_data$Sample) != benchmark_sample]

    for(samp in compared_samples){
      compdf <- .self$expression_data %>%
        filter(Sample == samp) %>% mutate(Exp2 = Exp) %>% select(!Sample & !Exp) %>%
        merge(d1, by = c('Symbol', 'Gene', 'Replicate')) %>%
        group_by(Symbol) %>%
        summarise(P = ifelse(n() > 1, t.test(Exp1, Exp2, paired = F, alternative = 'two.sided', var.equal = T)$p.value, NA),
                   Log2FC = log2(mean(Exp2)/mean(Exp1))) %>%
        mutate(Sample = samp)

      log2fc_pval_dataframe <<- rbind(.self$log2fc_pval_dataframe, compdf)

      .self$volcano_plots[[toString(samp)]] = ggplot(data = compdf, aes(x = Log2FC, y = -log10(P))) +
        geom_point() +
        theme_classic() +
        labs(title = samp)
    }

    log2fc_pval_dataframe <<- na.omit(log2fc_pval_dataframe)
  },

  GeneRegulation = function(log2fc_threshold = 1.5, pval_threshold = 0.01){
    goi_df <- .self$log2fc_pval_dataframe %>%
      filter((Log2FC < -abs(log2fc_threshold) | Log2FC > abs(log2fc_threshold)) & P < pval_threshold) %>%
      mutate(Direction = ifelse(Log2FC > 0, 1, 0))

    dir_counts <- goi_df %>%
      group_by(Sample) %>%
      summarise(Up = sum(Direction), Down = -(n() - sum(Direction))) %>%
      pivot_longer(!Sample, names_to = 'Direction', values_to = 'Sum')

    regulation_barplot <<- list('plot' = ggplot(dir_counts, aes(as.factor(Sample), Sum, fill = Direction)) +
      geom_bar(stat = 'identity', position = 'identity') +
      xlab("Sample")+ylab("Number of genes") +
      scale_fill_discrete(name = "",labels = c("Upregulated","Downregulated")) +
      labs(fill = "") +
      theme_classic() +
      theme(legend.position = "top"))

    goi_expression <<- .self$expression_data %>%
      filter(Symbol %in% goi_df$Symbol) %>%
      mutate(Sample = paste(Sample, Replicate, sep = '.')) %>%
      select(!Gene & !Replicate) %>%
      pivot_wider(values_from = Exp, names_from = Symbol) %>%
      column_to_rownames("Sample") %>%
      as.matrix()
  },

  PCA = function(clusters){

    elbow_plot <- fviz_nbclust(.self$goi_expression, kmeans, method = "wss") +
        labs(subtitle = "Elbow method")

    #ADD SILHOUETTE FUNCTION TO CALCULATE # CLUSTERS
    pca_plot <- autoplot(kmeans(.self$goi_expression, clusters),
              data = .self$goi_expression,
              label = TRUE,
              frame = TRUE,
              label.size = 3,
              colour = "cluster",
              label.repel= TRUE) +
       theme_classic() +
       theme(legend.position = "top")

    as_tibble(.self$goi_expression, rownames = "Sample")
    sample_pca<-prcomp(.self$goi_expression)
    pc_loadings<-sample_pca$rotation %>%
       as_tibble(rownames = "Gene")

    loadings_plot <- autoplot(ggplot(data = pc_loadings) +
       geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                    arrow = arrow(length = unit(0.1, "in")),
                    colour = "red") +
       geom_label_repel(aes(x = PC1, y = PC2, label = Gene)) +
       scale_x_continuous(expand = c(0.02, 0.02)) +
       labs(x = 'PC1', y = 'PC2') +
       theme_classic(), label.repel = T)

    pca <<- list(
        'pca_plot' = pca_plot,
       'elbow_plot' = elbow_plot,
       'loadings' = pc_loadings,
       'loadings_plot' = loadings_plot
    )
  },

  Pathways = function(method = 'both', top_path_n = 5){

    .self$pathways[['plots']] = list()

    path_barplot_wrapper <- function(pathdb, col){
      pathdb <- pathdb[order(pathdb[,col], decreasing = TRUE),][1:top_path_n,]

      length(pull(pathdb, col))

      ggplot(mapping = aes(x = reorder(pull(pathdb, 'Pathway'), -pull(pathdb, col)), y = pull(pathdb, col))) +
        geom_bar(stat = 'identity') +
        scale_y_continuous(expand = c(0, 0)) +
        xlab('Pathway') + ylab('Number of genes') +
        theme_classic() +
        theme(
          axis.text.x = element_text(size=12, angle=60,hjust=0.5,vjust=0.5, color = "black"),
          axis.ticks.x  = element_blank()
        )
    }

    ensembl_symbol_conversion <- .self$expression_data %>%
      select(Gene | Symbol) %>%
      unique()

    goi_long <- data.frame(.self$goi_expression) %>%
      rownames_to_column('Sample') %>%
      pivot_longer(!Sample, names_to = 'Symbol', values_to = 'Exp')

    if(method == 'both' | method == 'MitoCarta'){
      freq_temp <- read_sheet(.self$mitocarta_data$path, .self$mitocarta_data$pathway.sheet) %>%
        select(MitoPathway | Genes) %>%
        separate_rows(Genes, sep = ',') %>%
        na.omit() %>% unique() %>%
        mutate(Genes = gsub(" ", "", Genes)) %>%
        filter(Genes %in% goi_long$Symbol) %>%
        merge(.self$log2fc_pval_dataframe, by.x = 'Genes', by.y = 'Symbol') %>%
        mutate(Sig.Up = ifelse(Log2FC > 0, T, F)) %>%
        mutate(Sig.Down = ifelse(Log2FC < 0, T, F))

      freq2plot <- freq_temp %>%
        group_by(MitoPathway) %>%
        summarise(Sig.Up = sum(Sig.Up), Sig.Down = sum(Sig.Down)) %>%
        rename(Pathway = MitoPathway)

      .self$pathways[['plots']][['mitocarta.up']] = path_barplot_wrapper(freq2plot, 'Sig.Up')
      .self$pathways[['plots']][['mitocarta.down']] = path_barplot_wrapper(freq2plot, 'Sig.Down')

      .self$pathways[['mitocarta_frequencies']] = freq_temp %>%
        group_by(MitoPathway, Sample) %>%
        summarise(Sig.Up = sum(Sig.Up), Sig.Down = sum(Sig.Down))

      n_total_up = sum(freq_temp$Sig.Up)
      n_total_down = sum(freq_temp$Sig.Down)

      .self$pathways[['mitocarta_fisher_tests']] = freq_temp %>%
        group_by(MitoPathway) %>%
        summarise(P = fisher.test(
          matrix(
            c(sum(Sig.Up), sum(Sig.Down), n_total_up - sum(Sig.Up), n_total_down - sum(Sig.Down)),
            ncol = 2
          ),
          alternative = "two.sided"
        )$p.value)

      .self$pathways[['plots']][['heatmaps']] = list()
      top_paths <- append(freq2plot[order(freq2plot[,'Sig.Up'], decreasing = TRUE),][1:top_path_n,]$Pathway, freq2plot[order(freq2plot[,'Sig.Down'], decreasing = TRUE),][1:top_path_n,]$Pathway)
      top_pathway_genes_l2fcp <<- freq_temp %>%
        filter(Genes %in% goi_long$Symbol & MitoPathway %in% top_paths)

      hm_db <- .self$top_pathway_genes_l2fcp %>%
        select(Genes | Log2FC | MitoPathway | Sample) %>%
        mutate(Log2FC = as.numeric(Log2FC)) %>%
        pivot_wider(names_from = Sample, values_from = Log2FC)

      for(path in unique(top_pathway_genes_l2fcp$MitoPathway)){
        path_hmdb <- hm_db %>%
          filter(MitoPathway == path) %>%
          select(!MitoPathway) %>%
          column_to_rownames(var="Genes") %>%
          as.matrix()

        test_hmdb <<- path_hmdb

        path_hmdb <- path_hmdb[, unique(.self$log2fc_pval_dataframe$Sample)]

        .self$pathways[['plots']][['heatmaps']][[path]] = heatmaply(path_hmdb, Colv = NA, col = colorRampPalette(brewer.pal(3, "Blues"))(25))
      }

    }

    if(method == 'both' | method == 'Panther'){
      goi_logp <- .self$log2fc_pval_dataframe %>%
        merge(ensembl_symbol_conversion, by = 'Symbol') %>%
        filter(Symbol %in% goi_long$Symbol)

      .self$pathways[['panther_frequencies']] = .self$panther_query(filter(goi_logp, Log2FC > 0)$Gene) %>%
        mutate(Dir = 'Sig.Up') %>%
        rbind(
          .self$panther_query(filter(goi_logp, Log2FC < 0)$Gene) %>%
            mutate(Dir = 'Sig.Down')
        ) %>%
        pivot_wider(names_from = 'Dir', values_from = 'Fold.Enrichment')

      .self$pathways[['plots']][['panther.up']] = path_barplot_wrapper(.self$pathways[['panther_frequencies']], 'Sig.Up')
      .self$pathways[['plots']][['panther.down']] = path_barplot_wrapper(.self$pathways[['panther_frequencies']], 'Sig.Down')
    }
  }
  )
)
