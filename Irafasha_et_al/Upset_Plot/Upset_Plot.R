pacman::p_load(readxl, tidyverse, UpSetR, ComplexUpset, patchwork)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
# specifying the path for file
path <- "./"

# set the working directory 
setwd(path)

# accessing all the sheets 
sheet = excel_sheets("DESeq_DGE_results_Irafasha_et_al.xlsx")

# applying sheet names to dataframe names
dge_results = lapply(setNames(sheet, sheet), 
                    function(x){
                      read_excel("DESeq_DGE_results_Irafasha_et_al.xlsx", sheet=x)
                      })

# attaching all dataframes together
dge_results.df = bind_rows(dge_results, .id="treatment") %>% # bind sheets with sheet-name as ID under column "treatment"
  separate(treatment, c('treatment', 'time'), sep = "_") %>% # split treatment using underscore delimiter 
  as_tibble()

# printing data of all sheets
(dge_results.df)

## Identify treatments and timepoints
unique(dge_results.df$treatment)
unique(dge_results.df$time)

dge_results.signif <- dge_results.df %>% 
  filter(
    padj <= 0.05,
    abs(log2FoldChange) >= 2
         ) # filter adj p val <= .05 and absolute |logFC| >= 2

## Subset by timepoint 
dat.6hrs <- dge_results.signif %>%
  filter(time == "06HRS") %>% ## change for all the other time-points
  dplyr::select(1,3) %>% 
  reshape2::dcast(.,gene ~ treatment, fun.aggregate = function(x) 1L, fill = 0L) 

dat.12hrs <- dge_results.signif %>%
  filter(time == "12HRS") %>% ## change for all the other time-points
  dplyr::select(1,3) %>% 
  reshape2::dcast(.,gene ~ treatment, fun.aggregate = function(x) 1L, fill = 0L)

dat.18hrs <- dge_results.signif %>%
  filter(time == "18HRS") %>% ## change for all the other time-points
  dplyr::select(1,3) %>% 
  reshape2::dcast(.,gene ~ treatment, fun.aggregate = function(x) 1L, fill = 0L)

## Library upsetR Try for 6h
UpSetR::upset(dat.6hrs,
      sets = c("GR24", "IS2730", "SRN39", "IS41724", "IS27146"),
      # sets = set_vars,
      mainbar.y.label = "Number of genes intersections", 
      sets.x.label = "Counts by Condition",
      order.by = "freq",
      )

## A fancy update on the plots with library complex-upset
p1 = upset(
  dat.6hrs,
  c('GR24', 'IS2730', 'SRN39', 'IS41724', 'IS27146'),
  queries=list(
    upset_query(set = 'GR24', fill='orange'),
    upset_query(set = 'IS2730', fill='blue'),
    upset_query(set = 'SRN39', fill='coral3'),
    upset_query(set = 'IS41724', fill='navy'),
    upset_query(set = 'IS27146', fill='black')
  ),
  base_annotations=list(
    '#gene intersections'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      # +theme_classic()
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colored
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    ) 
    + ylab("# DEGs at 06 HPI")

  ),
  sort_sets='ascending',
  sort_intersections='descending'
)

## 12 h
p2 = upset(
  dat.12hrs,
  c('GR24', 'IS2730', 'SRN39', 'IS41724', 'IS27146'),
  queries=list(
    upset_query(set = 'GR24', fill='orange'),
    upset_query(set = 'IS2730', fill='blue'),
    upset_query(set = 'SRN39', fill='coral3'),
    upset_query(set = 'IS41724', fill='navy'),
    upset_query(set = 'IS27146', fill='black')
  ),
  base_annotations=list(
    '#gene intersections'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      # +theme_classic()
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colored
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    ) 
    + ylab("# DEGs at 12 HPI")
    
  ),
  sort_sets='ascending',
  sort_intersections='descending'
)

## 18 h
p3 = upset(
  dat.18hrs,
  c('GR24', 'IS2730', 'SRN39', 'IS41724', 'IS27146'),
  queries=list(
    upset_query(set = 'GR24', fill='orange'),
    upset_query(set = 'IS2730', fill='blue'),
    upset_query(set = 'SRN39', fill='coral3'),
    upset_query(set = 'IS41724', fill='navy'),
    upset_query(set = 'IS27146', fill='black')
  ),
  base_annotations=list(
    '#gene intersections'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      # +theme_classic()
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colored
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    ) 
    + ylab("# DEGs at 18 HPI")
    
  ),
  sort_sets='ascending',
  sort_intersections='descending'
)

p <- (p1/p2/p3)

ggsave(plot = p, filename = "upsetRplots.svg", height = 20, width = 18) ## play around with height/width to get the desired plot dimensions

