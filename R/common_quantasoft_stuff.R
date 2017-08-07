Quantaimport <- function(){
  require(dplyr)
  require(readr)
  filenames <- list.files(pattern="*.csv")
  bigfile <- temp<- data.frame()

  #Just select desired columns rather than figure out how/why excel shat out an extra comma
  for(file in filenames){
    temp <- read_csv(file)
    #Pick a pretty specific colname to check if it's a quantasoft export
    if('PoissonConfMin' %in% colnames(temp)){
      temp <- select(temp, Well:PoissonFractionalAbundanceMin68)
      temp$Plate = file
      # #If you want to add anything cute from the filename, now is the time to do it, i.e
      # # temp$DR <- strsplit(file, '_')[[1]][6]
      # # temp$User <- strsplit(file, '_')[[1]][4]
      # # temp$Day <- strsplit(file, '_')[[1]][5]
      #
      # if(grepl('Kiyoung', file)){
      #   temp$User <- 'Kiyoung'
      # }
      # if(grepl('Prasanthi', file)){
      #   temp$User <- 'Prasanthi'
      # }
      # if(grepl('Nathan', file)){
      #   temp$User <- 'Nathan'
      # }
      # if(grepl('Tiffany', file)){
      #   temp$User <- 'Tiffany'
      # }
      bigfile <- rbind(bigfile, temp)
      cat(c('finished reading in', file, "\n"))
    } else {
      cat(c('skipping', file, "\n"))
    }
  }
  return(bigfile)
}
#Grab all wells from a plate with only singles. Pick 2 adjacent wells and merge them.
#Can also shuffle tests if you want that.
merge_plate <- function(myplate, n_to_combine = 2, shuffle = FALSE){
  require(dplyr)
  if('M' %in% substr(myplate$Well, 1, 1)){
   stop('NO MERGED WELLS ALLOWED GO AWAY')
  }

  temp <- ''
  samplelist <- unique(myplate$Sample)
  targetlist <- unique(myplate$Target)
  #mergedwells is output, merged counter keep track of where to put next row
  #Done this way because I think rbind is slow? I dunno.
  mergedwells <- myplate[1,]
  #merged_counter does something important maybe.
  merged_counter <- 1

  for(target in targetlist){
    targetplate <- filter(myplate, Target == target)

    for(sample in samplelist){
      #Exact match for sample only
      temp <- paste('^', sample, '$', sep = '')
      well_list <- grep(temp, targetplate$Sample)
      well_counter <- 1
      if(length(well_list) %% n_to_combine == 1){
        stop('Plate does not divide evenly with this n_to_combine')
        }

      #shuffle order at random
      if(shuffle == TRUE){
        well_list <- sample(well_list)
      }
      #loops through sample list, should stop when well counter == length of well list
      while(well_counter < length(well_list)){
        tomerge <- data.frame()
        for(x in 1:n_to_combine -1){
          #Makes a new data frame of target wells to combine. Probably not efficient.
          #This breaks if the well list doesn't divide evenly Stick and IF
          if(well_counter + x <= length(well_list)){
            tomerge <- rbind(tomerge, targetplate[well_list[well_counter +x],])
          }
        } #end of n to combine loop
        #Recruit another function because this one is getting messy.
        mergedwells[merged_counter,] <- merge_wells(tomerge)
        #Merged counter tells us where we are in the output data frame
        #Well counter tells us where we are in the list of wells.
        merged_counter <- merged_counter + 1
        well_counter <- well_counter + n_to_combine
      } #end of well list loop
    } #end of sample loop
  } #end of target loop
  return(mergedwells)
} #end of function
#combine 2 wells. Lazy.

merge_wells <- function(tomerge){
  #Set colnames and size of metawell in horrible wrong way
  metawell <- tomerge[1,]
  metawell[1,] <- NA

  #These 3 are really all that matters.
  metawell$AcceptedDroplets <- sum(tomerge[,'AcceptedDroplets'])
  metawell$Positives <- sum(tomerge[,'Positives'])
  metawell$Negatives <- sum(tomerge[,'Negatives'])
  #I expect this will only ever be one value for these, but make sure
  #Don't really care about most of these, but why not.
  metawell$Well <- paste(tomerge$Well, collapse = ',')
  metawell$ExptType <- paste(unique(tomerge$ExptType), collapse = ',')
  metawell$Experiment <- paste(unique(tomerge$Experiment), collapse = ',')
  metawell$Sample <- paste(unique(tomerge$Sample), collapse = ',')
  metawell$TargetType <- paste(unique(tomerge$TargetType), collapse = ',')
  metawell$Target <- paste(unique(tomerge$Target), collapse = ',')
  metawell$Status <- paste(unique(tomerge$Status), collapse = ',')
  #Droplet size is 0.00085.
  #Quantasoft rounds to 2 decimals
  #CpD = -ln(negatives/total)
  #Conc = CpD / 0.00085
  metawell$Concentration <- round(-log(metawell$Negatives / metawell$AcceptedDroplets) / 0.00085, 2)
  metawell$Supermix <- paste(unique(tomerge$Supermix), collapse = ',')
  metawell$CopiesPer20uLWell <- metawell$Concentration * 20
  # metawell$DR <- paste(unique(tomerge$DR), collapse = ',')
  # metawell$User <- paste(unique(tomerge$User), collapse = ',')
  # metawell$Day <- paste(unique(tomerge$Day), collapse = ',')
  #Gotta check whether we're ch1 or ch2 for ratio
  if(metawell$TargetType == 'Ch1Unknown'){
    otherneg <- sum(as.numeric(tomerge$`Ch1+Ch2-` + tomerge$`Ch1-Ch2-`))
    otherconc <- round(-log(otherneg / metawell$AcceptedDroplets) / 0.00085, 2)
    metawell$Ratio <- signif(metawell$Concentration / otherconc, 3)
  }

  if(metawell$TargetType == 'Ch2Unknown'){
    otherneg <- sum(as.numeric(tomerge$`Ch1-Ch2+` + tomerge$`Ch1-Ch2-`))
    otherconc <- round(-log(otherneg / metawell$AcceptedDroplets) / 0.00085, 2)
    metawell$Ratio <- signif(otherconc / metawell$Concentration, 3)
  }
  #And then some other stuff I don't care about
  return(metawell)
}

heatmap <- function(bigfile){
  require(ggplot2)
  require(dplyr)
  #Pull Row, column data
  bigfile$Row <- substr(bigfile$Well, 1,1)
  bigfile$Column <-substr(bigfile$Well, 2,3)
  #No merged wells allowed, Only take hex channel to avoid duplicate wells
  singles <- filter(bigfile, Row != 'M', TargetType =='Ch2Unknown')
  #Geom tile makes heat maps, hooray.
  #Can modify color with scale_fill_discrete() or whatever
  #Fills in blanks on plots
   heatmap <- ggplot(singles, aes(x = Column, y = Row)) +
    geom_tile(aes(fill = AcceptedDroplets)) +
    scale_x_discrete(position = 'top',
                     limits = c('01','02','03', '04', '05', '06',
                                '07', '08', '09', '10', '11', '12')) +
    scale_y_discrete(limits = c('H','G','F','E','D','C','B','A')) +
    scale_fill_gradient2(limits = c(0,30000), midpoint = 18000) +
    labs(title = 'Droplet counts', x = 'Column', y = 'Row')

  cat('Median Drops', median(bigfile$AcceptedDroplets),
      '\n<10k', sum(singles$AcceptedDroplets < 10000),
      '\n<12k', sum(singles$AcceptedDroplets < 12000))
  return(heatmap)
}
failmap <-function(bigfile){
  dropsum <- filter(bigfile, TargetType == 'Ch2Unknown') %>%
    group_by(Well) %>%
    summarise(drops = mean(AcceptedDroplets),
              fail10k = sum(AcceptedDroplets <10000),
              fail10krate = fail10k/n(),
              count = n())

  dropsum$Row <- substr(dropsum$Well, 1,1)
  dropsum$Column <-substr(dropsum$Well, 2,3)

  failmap <- ggplot(data = dropsum) + geom_tile(aes(x = Column, y = Row, fill = fail10krate)) +
    scale_x_discrete(position = 'top',
                     limits = c('01','02','03', '04', '05', '06',
                                '07', '08', '09', '10', '11', '12')) +
    scale_y_discrete(limits = c('H','G','F','E','D','C','B','A'))
  return(dropsum)
}

#Just calls the summary function with a bunch of stats. Saves some typing
imlazy <- function(plate){
  require(dplyr)
  ABLsummary<-filter(plate, TargetType == 'Ch2Unknown') %>%
    group_by(Sample) %>%
     summarise(ratio_mean = mean(Ratio),
              PercentIS = (ratio_mean * 0.93 * 100),
              MR = -log10(PercentIS/100),
              ratio_cv = (sd(Ratio) / mean(Ratio) * 100),
              ABL = mean(Concentration),
              ABL_CV = (sd(Concentration) / mean(Concentration) *100),
              Events = mean(AcceptedDroplets),
              ABLcopies = as.integer(ABL * Events * 0.00085),
              N_tests = n())

  return(ABLsummary)
}

#Somebody else's function to put linear models on ggplots
stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = 'identity',
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}

library(ggplot2)
StatSmoothFunc <- ggproto("StatSmooth", Stat,

                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))

                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }

                            params
                          },

                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }

                            if (is.null(data$weight)) data$weight <- 1

                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }

                            if (is.character(method)) method <- match.fun(method)

                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))

                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                             list(a = format(coef(m)[1], digits = 3),
                                                  b = format(coef(m)[2], digits = 3),
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))

                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)

                          },

                          required_aes = c("x", "y")
)

shitty_well_filter <- function(data){
  require(dplyr)
  return(filter(data, AcceptedEvents > 10000))
}

replace_missing_merged <- function(d){
  #Pull off the merged and singles into different sets.
  merged <- d[grep('M', d$Well), ]
  singlets <- d[grep("[A-H]", d$Well), ]

  #Adds rows that exist in singlets, but were not included in merged wells
  if(nrow(singlets) == 0){
    print (c(file,'Has no single well data, here, you can have it back'))
    return(d)
  }
  if(nrow(merged) == 0){
    print( c(file, 'Has no merged data. Did you mean to use merge_plate?'))
    return(d)
  }

  for(x in (1:nrow(singlets))){
    #Reads as "If row x well entry is not in mergedwells, return True
    temp <- singlets[x,'Well']
    if( !(TRUE %in% grepl(temp, merged$MergedWells))){
      print (c(singlets$Sample[x], 'Not in merged dataset, adding'))
      merged <- rbind(merged, singlets[x,])
    }
  }
  return(merged)
}

Use_lm_function <- function(model, data){
  terms <- data.frame(data)
  terms[,2] <- coef(model)[1]
  for(x in 2:length(coef(model))){
    terms[,x+1] <- (coef(model)[x] * (data^(x-1)))
  }
  terms[,'data'] <- NULL
  return(rowSums(terms))
}

Alternate_ratio_calc <- function(data){
  require(dplyr)
  #Filter out ABLS, add their concentration to BCRABLS
  ABLs <- filter(data, TargetType =='Ch2Unknown')
  WithABL <- merge(x = data, y = select(ABLs, Well, ABLCPM = Concentration), by.x = 'Well', by.y = 'Well')
  BCRABL <- filter(WithABL, TargetType =='Ch1Unknown')
  #Add new ratio to BCRABLs, extend to whole dataset new ratio = BCRABL/BCRABL+ABL
  BCRABL$Ratiofix <- BCRABL$Concentration/ (BCRABL$ABLCPM -BCRABL$Concentration)
  data <- merge(x = data, y = select(BCRABL, Well, Ratiofix), by.x = 'Well', by.y = 'Well')
  #Column switcheroo
  data$oldRatio <- data$Ratio
  data$Ratio <- data$Ratiofix
  data$Ratiofix <- NULL
  return(data)
}


