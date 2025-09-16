#############-- Install packages --#############
# install.packages("readxl")
# install.packages("ggplot2")
# install.packages("ggthemes")
# install.packages("writexl")
# install.packages("stringr")
# install.packages("tidyverse")
# BiocManager::install("ComplexHeatmap")
# install.packages("circlize")

#############-- Libraries --#############
library(readxl)
library(ggplot2)
library(writexl)
library(stringr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

#############-- Functions --#############

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

theme_Publication <- function(base_size = 14) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

read_in_raw_bli_data <- function(file, results_folder) {
    total_antibodies <- c()
    sheet_names <- excel_sheets(file)
    list_of_sheets <- lapply(sheet_names, function(x) read_excel(file, sheet = x, skip=1))
    names(list_of_sheets) <- sheet_names
    for (i in names(list_of_sheets)) {
        if (!(i %in% c("Run9", "Run10", "Run11"))) {
            df <-  as.data.frame(list_of_sheets[[i]])
            # if ("1A8" %in% colnames(df)) {
            #     print(i)
            # }
            df <- df[!(df$Sample %in% c("1a8", "1A8", "1A8 -new", "1A8 - new","1A8-NEW","66", "64")), !(colnames(df) %in% c("1a8", "1A8","1A8 -new", "1A8 - new","1A8-NEW", "66", "64")) ]
            list_of_sheets[[i]] <- df
            list_of_sheets[[i]] <- list_of_sheets[[i]][list_of_sheets[[i]]$Sample!="no mAb",]
            row.names(list_of_sheets[[i]]) <- list_of_sheets[[i]]$Sample
            list_of_sheets[[i]]$Sample <- NULL
        }  else { 
            list_of_sheets[[i]] <- NULL
        }
    }
    norm_list <- list_of_sheets
    for (i in names(norm_list)) {      
        print(i)
        norm_list[[i]] <- as.matrix(as.data.frame(norm_list[[i]]))
        class(norm_list[[i]]) <- "numeric"
        norm_antibody <- colnames(norm_list[[i]] )
        norm_antibody <- norm_antibody[!(norm_antibody %in% c("Sample", "no mAb" ))]
        for (a in norm_antibody) {
            norm_list[[i]][, a] <- norm_list[[i]][, a]/norm_list[[i]][, "no mAb"] 
        }
        norm_list[[i]] <- as.data.frame(norm_list[[i]][, !(colnames( norm_list[[i]]) %in% c("no mAb"))])
        norm_list[[i]]$Sample <- toupper(rownames(norm_list[[i]]))
        colnames(norm_list[[i]]) <- toupper(colnames(norm_list[[i]]) )
        colnames(norm_list[[i]]) <- str_remove_all(colnames(norm_list[[i]]), " ")
        norm_list[[i]]$SAMPLE <-  str_remove_all(norm_list[[i]]$SAMPLE, " ")
        colnames(norm_list[[i]]) <- str_remove_all(colnames(norm_list[[i]]), "-NEW2")
        norm_list[[i]]$SAMPLE <-  str_remove_all(norm_list[[i]]$SAMPLE, "-NEW2")
        norm_antibody <- colnames(norm_list[[i]] )
        total_antibodies <- c(total_antibodies,norm_antibody)
    }
    write_xlsx(norm_list, file.path(results_folder, "mABNormdata.xlsx"))
    return(list("values"=norm_list, "total_antibodies"=total_antibodies))
}

generate_pair_combinations <- function(orderp, orders) {
    dfcomb <- as.data.frame(matrix(nrow=(length(primary)*(length(secondary)-1)), ncol=3))
    colnames(dfcomb) <- c("combpair", "primary", "secondary")
    counter=1
    counterpair=1
    tmp <- c()
    for (p in orderp) { 
        for (s in orders) { 
            if (p!=s) {
                if (!(paste0(p,"-",s) %in% tmp)){ 
                    dfcomb[counter,"combpair"] <- paste0("comb", counterpair)
                    dfcomb[counter,"primary"] <- p
                    dfcomb[counter,"secondary"] <- s
                    counter <- counter+1
                    dfcomb[counter,"combpair"] <- paste0("comb", counterpair)
                    dfcomb[counter,"primary"] <- s
                    dfcomb[counter,"secondary"] <- p
                    counter <- counter+1
                    counterpair <- counterpair+1
                    tmp <- c(tmp, paste0(p,"-",s))
                    tmp <- c(tmp, paste0(s,"-",p))
                }
            }
        }
    }
    return(dfcomb)
}

organize_data <- function(data,dfcomb, primary, secondary) {
    ### Generate averages across all normalized data for each pair 
    dfpairsallnormvalues <- list()
    datapairs <- list()
    datapairsNoNegVals <- list()
    datapairsZeros <- list()

    for (p in primary) {
        tmp <- c()
        for (s in secondary){
            for (i in names(data$values)) {
                if (((p %in% colnames(data$values[[i]])) & (s %in% data$values[[i]]$SAMPLE))) {
                    tmp <- rbind(tmp, c(s, data$values[[i]][data$values[[i]]$SAMPLE==s, p]))
                }
            }
        }
        tmp <- as.data.frame(tmp)
        colnames(tmp) <- c("secondary", "values")
        tmp$values <- as.numeric(tmp$values)
        dfpairsallnormvalues[[p]] <- tmp

        tmpNoNeg <- tmp[tmp$values >0, ]
        tmpZeros <- tmp
        tmpZeros[tmpZeros < 0] <- 0

        tmpavg <- tmp %>% group_by(secondary) %>% summarise(mean = mean(values))
        tmpavg <- as.data.frame(tmpavg)
        datapairs[[p]] <- tmpavg

        tmpavg <- tmp %>% group_by(secondary) %>% summarise(mean = mean(values))
        tmpavg <- as.data.frame(tmpavg)
        datapairs[[p]] <- tmpavg

        tmpavg <- tmpNoNeg %>% group_by(secondary) %>% summarise(mean = mean(values))
        tmpavg <- as.data.frame(tmpavg)
        datapairsNoNegVals[[p]] <- tmpavg

        tmpavg <- tmpZeros %>% group_by(secondary) %>% summarise(mean = mean(values))
        tmpavg <- as.data.frame(tmpavg)
        datapairsZeros[[p]] <- tmpavg
    }

    # box plot of all values 
    dfpairsallnormvaluesdf <- c()
    for (i in names(dfpairsallnormvalues)) {
        dfpairsallnormvalues[[i]]$primary <- rep(i, nrow(dfpairsallnormvalues[[i]]))
        dfpairsallnormvaluesdf <- rbind(dfpairsallnormvaluesdf, dfpairsallnormvalues[[i]] )
    }
    pair <- c()
    combid <- c()
    for (i in rownames(dfpairsallnormvaluesdf)) {
        if (dfpairsallnormvaluesdf[i, "secondary"]==dfpairsallnormvaluesdf[i, "primary"]) {
            pair <- c(pair, "exactpair")
            combid <- c(combid, paste0("comb", length(unique(dfcomb$combpair))+1))
        } else { 
            pair <- c(pair, paste(dfpairsallnormvaluesdf[i, "primary"], dfpairsallnormvaluesdf[i, "secondary"], sep="-"))
            combid <- c(combid, dfcomb[(dfcomb$primary==dfpairsallnormvaluesdf[i, "primary"] & dfcomb$secondary==dfpairsallnormvaluesdf[i, "secondary"]), "combpair"]) 
        }
    }
    dfpairsallnormvaluesdf$pair <- pair
    dfpairsallnormvaluesdf$combid <- combid
    dfpairsallnormvaluesdf$combid <- factor(dfpairsallnormvaluesdf$combid, 
        levels=c(unique(dfcomb$combpair), paste0("comb", length(unique(dfcomb$combpair))+1)))
    ggplot(dfpairsallnormvaluesdf, aes(x=combid, y=values, group=pair)) +
        geom_boxplot(position=position_dodge(1)) + theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5), 
            axis.title.x=element_blank()) + labs(y="norm values")
    ggsave(file.path(results_folder, "boxplot_normvalues.png"), width=8, height=4, dpi=300)

    ggplot(dfpairsallnormvaluesdf[dfpairsallnormvaluesdf$values >0, ], aes(x=combid, y=values, group=pair)) +
        geom_boxplot(position=position_dodge(1)) + theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
        axis.title.x=element_blank()) + labs(y="norm values")
    ggsave(file.path(results_folder, "boxplot_normvaluesNoZeros.png"), width=8, height=4, dpi=300)

    ggplot(dfpairsallnormvaluesdf, aes(x=values)) + 
    geom_density() + theme_Publication() + labs(x="norm values")
    ggsave(file.path(results_folder, "density_normvalues.png"), width=4, height=4, dpi=300)

    ggplot(dfpairsallnormvaluesdf[dfpairsallnormvaluesdf$values>0, ], aes(x=values)) + 
    geom_density() + theme_Publication() + labs(x="norm values")
    ggsave(file.path(results_folder, "density_normvalues_noZeros.png"), width=4, height=4, dpi=300)
    return(list("allnormvalues"=dfpairsallnormvaluesdf,"normavgs"=datapairs,
        "normavgsNoNegValues"=datapairsNoNegVals, "normavgsZeros"=datapairsZeros))
}

absolute_diff_between_pairs <- function(df,label) { 
    dfmelt <- c()
    for (i in names(df)) {
        tmp <- df[[i]]
        tmp$primary <- rep(i, nrow(tmp))
        dfmelt <- rbind(dfmelt, tmp)
    }
    dfmelt$mean <- dfmelt$mean*100
    dfabspairs <- as.data.frame(matrix(nrow=length(unique(dfcomb$combpair)), ncol=4))
    colnames(dfabspairs) <- c("combpair", "pair", "absdiff", "avg")
    counter=1
    for (v in unique(dfcomb$combpair)) {
        p= dfcomb[dfcomb$combpair==v, "primary"]
        if ("EM4C04" %in% p) {
            absdiff<-0
            test <- dfmelt[dfmelt$primary=="EM4C04" ,]
            if (p[2] %in% test$secondary) {
                avg <- dfmelt[dfmelt$primary=="EM4C04" & dfmelt$secondary==p[2], "mean"]
           } else {
            avg=NA
        }
        } else {
            absdiff <- abs(dfmelt[dfmelt$primary==p[1] & dfmelt$secondary==p[2], "mean"] - dfmelt[dfmelt$primary==p[2] & dfmelt$secondary==p[1], "mean"])
            avg <- mean(c(dfmelt[dfmelt$primary==p[1] & dfmelt$secondary==p[2], "mean"], dfmelt[dfmelt$primary==p[2] & dfmelt$secondary==p[1], "mean"]))
        }
        dfabspairs[counter, "combpair"] <- v
        dfabspairs[counter, "pair"] <- paste(p[1],p[2],sep="-")
        if(identical(absdiff, numeric(0))) {
            dfabspairs[counter,"absdiff"] <- NA
            dfabspairs[counter,"avg"] <- NA
        } else { 
            dfabspairs[counter,"absdiff"] <- absdiff
            dfabspairs[counter,"avg"] <- avg 
        }
        counter = counter+1
    }

    dfabspairs <- dfabspairs[rev(order(dfabspairs$absdiff)), ]
    dfabspairstmp <- dfabspairs[complete.cases(dfabspairs), ]
    dfabspairstmp$pair <- factor(dfabspairstmp$pair, levels=dfabspairstmp$pair)
    ggplot(dfabspairstmp, aes(x=pair, y=absdiff)) +
        geom_bar(stat="identity") + theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5)) +
        labs(y="abs difference primary vs secondary &\nsecondary vs primary",x ="pairs")
    ggsave(file.path(results_folder, paste0("absolute_pair_diff_",label,".png")), width=8,height=4, dpi=500)

    return(dfabspairs)
}

generate_initial_avgs_hm <- function(datapairs, primary, secondary, label) {
    df <- as.data.frame(matrix(nrow=length(secondary), ncol=length(primary)))
    colnames(df) <- primary
    rownames(df) <- secondary
    for (i in names(datapairs)) {
        for (s in datapairs[[i]]$secondary) {
            df[s,i] <- datapairs[[i]][datapairs[[i]]$secondary==s,"mean"]   
        } 

    }
    df <-df[orders,orderp]
    df <- df*100
    write.csv(df, file.path(results_folder, paste0("averagedNormPairs_", label,".csv")))
    df <- as.matrix(data.frame(df))
    rownames(df) <- colnames(df)[1:nrow(df)]
    colnames(df) <- str_replace_all(colnames(df), "X", "mAb-")
    rownames(df) <- str_replace_all(rownames(df), "X", "mAb-")

    col_fun <- colorRamp2(c(-20, 0, 20, 40, 60, 80), c("darkred","red", "orange", "yellow","green", "darkgreen"))
    Stagecol <- factor(rep(c("695", "684","1303","386","385",  "VRCO1", "EM4C04"), c(3, 4, 1, 2, 3, 1, 1)),
    levels =c("695", "684","1303","386","385",  "VRCO1", "EM4C04"))
    Stagerow <- factor(rep(c("695", "684","1303","386","385", "VRCO1"), c(3, 4,1, 2, 3, 1)),
    levels =c("695", "684","1303","386","385", "VRCO1"))
    hmC <- Heatmap(df,
        name = "% residual\nbinding",  na_col = "white",
        col = col_fun,
        show_row_names = TRUE,
        row_split = Stagerow,
        cluster_row_slices = T,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_split = Stagecol,
        cluster_column_slices = F,
        column_title_rot = 90,
        row_title_rot = 0,
        border = F, width = unit(12, "cm"),height = unit(10, "cm"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
        if (!is.na(df[i,j])) {
            grid.text(round(df[i, j],0), x, y, gp=gpar(fontsize=8))
            }
        },
        heatmap_legend_param = list(
            legend_height = unit(3, "cm"),
            labels_gp = gpar(fontsize = 14),
            title_gp = gpar(fontsize = 14, fontface = "bold")
        ), use_raster = TRUE,  raster_quality = 5
    )
    png(file.path(results_folder, paste0("averagedNormPairs", label, ".png")), width = 8, height = 8, units = "in", res = 700)
    # pdf(file.path(deresults_path, "ComplexHeatmap.pdf")) # , width = 6, height = 6, units = "in", res = 700)

    heatmapCOV <- draw(hmC,
        heatmap_legend_side = "left"
    )
    dev.off()
    return(df)
}

generate_final_avgs_hm <- function(dfall, primary, secondary, label) {
    dffinal <- as.data.frame(matrix(nrow=length(secondary), ncol=length(primary)))
    colnames(dffinal) <- primary
    rownames(dffinal) <- secondary
    colnames(dfall) <- str_remove_all(colnames(dfall), "mAb-")
    rownames(dfall) <- str_remove_all(rownames(dfall), "mAb-")
    dffinal <- dffinal[orders,orderp]

    for (i in unique(dfcomb$combpair)) {
        dftmp <- dfcomb[dfcomb$combpair==i,]
        if ("EM4C04" %in% dftmp[,"primary"]) {
            dffinal[dftmp[,"primary"][2] , "EM4C04"] <- dfall[dftmp[,"primary"][2], "EM4C04"]
        } else { 
            if(is.na(dfall[dftmp[1,"primary"] ,dftmp[1,"secondary"]]) & is.na(dfall[dftmp[2,"primary"] ,dftmp[2,"secondary"]])) {
                dffinal[dftmp[1,"primary"] ,dftmp[1,"secondary"]] <- NA
                dffinal[dftmp[2,"primary"] ,dftmp[2,"secondary"]] <- NA
            } else if (!(is.na(dfall[dftmp[1,"primary"] ,dftmp[1,"secondary"]])) & is.na(dfall[dftmp[2,"primary"] ,dftmp[2,"secondary"]])) {
                dffinal[dftmp[1,"primary"] ,dftmp[1,"secondary"]] <- dfall[dftmp[1,"primary"] ,dftmp[1,"secondary"]]
                dffinal[dftmp[2,"primary"] ,dftmp[2,"secondary"]] <- dfall[dftmp[1,"primary"] ,dftmp[1,"secondary"]]
            } else if (is.na(dfall[dftmp[1,"primary"] ,dftmp[1,"secondary"]]) & !(is.na(dfall[dftmp[2,"primary"] ,dftmp[2,"secondary"]]))) {
                dffinal[dftmp[1,"primary"] ,dftmp[1,"secondary"]] <- dfall[dftmp[2,"primary"] ,dftmp[2,"secondary"]]
                dffinal[dftmp[2,"primary"] ,dftmp[2,"secondary"]] <- dfall[dftmp[2,"primary"] ,dftmp[2,"secondary"]]     
            } else {
                totalavg <- mean(c(dfall[dftmp[1,"primary"] ,dftmp[1,"secondary"]], dfall[dftmp[2,"primary"] ,dftmp[2,"secondary"]]))
                dffinal[dftmp[1,"primary"] ,dftmp[1,"secondary"]] <- totalavg
                dffinal[dftmp[2,"primary"] ,dftmp[2,"secondary"]] <- totalavg
            }
        }
    }

    write.csv(dffinal, file.path(results_folder, paste0("FINALaveragedNormPairs_", label,".csv")))
    dffinal <- as.matrix(data.frame(dffinal))
    rownames(dffinal) <- colnames(dffinal)[1:nrow(dffinal)]
    colnames(dffinal) <- str_replace_all(colnames(dffinal), "X", "mAb-")
    rownames(dffinal) <- str_replace_all(rownames(dffinal), "X", "mAb-")
    v=1
    for (c in colnames(dffinal)) {
        if(length(rownames(dffinal)[v:length(rownames(dffinal))])<2) {
            break
        }
        for(r in rownames(dffinal)[v:length(rownames(dffinal))]) {
            dffinal[r,c] <- NA
        }
        v=v+1
    }
    col_fun <- colorRamp2(c(0, 20, 40, 60, 80, 100), c("darkred","red", "orange", "yellow","green", "darkgreen"))

    Stagecol <- factor(rep(c("695", "684","1303","386","385",  "VRCO1", "EM4C04"), c(3, 4, 1, 2, 3, 1, 1)),
    levels =c("695", "684","1303","386","385",  "VRCO1", "EM4C04"))
    Stagerow <- factor(rep(c("695", "684","1303","386","385", "VRCO1"), c(3, 4,1, 2, 3, 1)),
    levels =c("695", "684","1303","386","385", "VRCO1"))
    hmC <- Heatmap(dffinal,
        name = "% residual\nbinding", na_col = "white",
        col = col_fun, 
        show_row_names = TRUE,
        row_split = Stagerow,
        cluster_row_slices = T,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_split = Stagecol,
        cluster_column_slices = F,
        column_title_rot = 90,
        row_title_rot = 0,
        border = F, width = unit(12, "cm"),height = unit(10, "cm"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
        if (!is.na(dffinal[i,j])) {
            grid.text(round(dffinal[i, j],0), x, y, gp=gpar(fontsize=8))
            }
        },
        heatmap_legend_param = list(
            legend_height = unit(3, "cm"),
            labels_gp = gpar(fontsize = 14),
            at = c(0, 20, 40, 60, 80, 100),
            labels=c(0, 20, 40, 60, 80, 100),
            #labels = c("", "Strong Comp.", "", "Medium Comp.", "No Comp.", ""),
            title_gp = gpar(fontsize = 14, fontface = "bold")
        ), use_raster = TRUE,  raster_quality = 5
    )
    png(file.path(results_folder, paste0("FINALaveragedNormPairs", label, ".png")), width = 8, height = 8.5, units = "in", res = 700)
    heatmapCOV <- draw(hmC,
        heatmap_legend_side = "left"
    )
    dev.off()
    pdf(file.path(results_folder, paste0("FINALaveragedNormPairs", label, ".pdf")), width = 8, height = 8.5)
    heatmapCOV <- draw(hmC,
        heatmap_legend_side = "left"
    )
    dev.off()
    col_fun2 <- colorRamp2(c(0, 100), c("#895129", "#026464"))
    dffinalbinary <- dffinal
    dffinalbinary[dffinalbinary > 50] <- 100
    dffinalbinary[dffinalbinary <= 50] <- 0

    hmC <- Heatmap(dffinalbinary,
        name = "% residual\nbinding", na_col = "white",
        col = col_fun2, 
        show_row_names = TRUE,
        row_split = Stagerow,
        cluster_row_slices = T,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_split = Stagecol,
        cluster_column_slices = F,
        column_title_rot = 90,
        row_title_rot = 0,
        border = F, width = unit(12, "cm"),height = unit(10, "cm"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
        if (!is.na(dffinal[i,j])) {
            grid.text(round(dffinal[i, j],0), x, y, gp=gpar(fontsize=10))
            }
        },
        heatmap_legend_param = list(
            legend_height = unit(3, "cm"),
            labels_gp = gpar(fontsize = 14),
            at = c(0, 20, 40, 60, 80, 100),
            labels=c(0, 20, 40, 60, 80, 100),
            #labels = c("", "Strong Comp.", "", "Medium Comp.", "No Comp.", ""),
            title_gp = gpar(fontsize = 14, fontface = "bold")
        ), use_raster = TRUE,  raster_quality = 5
    )
    png(file.path(results_folder, paste0("FINALaveragedNormPairs", label, "_binary.png")), width = 8, height = 8.5, units = "in", res = 700)
    heatmapCOV <- draw(hmC,
        heatmap_legend_side = "left"
    )
    dev.off()
    pdf(file.path(results_folder, paste0("FINALaveragedNormPairs", label, "_binary.pdf")), width = 8, height = 8.5)
    heatmapCOV <- draw(hmC,
        heatmap_legend_side = "left"
    )
    dev.off()
}

generate_VRC01figs <- function(df, label) {
    dftmp <- df[complete.cases(df), ]
    dftmpstmpVRC01 <- dftmp %>%
    filter(str_detect(pair, "VRC01"))
    if (TRUE %in% str_detect(dftmpstmpVRC01$pair, "1A8")) {
        clonotype_info <- c("684", "684", "684", "684", "1303", "695", "695", "695","386", "386","385", "385", "385", NA)
        order_pairs <- c("18-VRC01","32-VRC01", "1A8-VRC01", "1G3-VRC01",
            "76-VRC01", "4-VRC01", "70-VRC01", "12-VRC01",
            "16-VRC01", "14-VRC01", "6-VRC01", "10-VRC01", "28-VRC01", "EM4C04-VRC01")
    } else {
           clonotype_info <- c("684", "684", "684", "1303", "695", "695", "695","386", "386","385", "385", "385",  NA)
            order_pairs <- c("18-VRC01","32-VRC01", "1G3-VRC01",
            "76-VRC01", "4-VRC01", "70-VRC01", "12-VRC01",
            "16-VRC01", "14-VRC01", "6-VRC01", "10-VRC01", "28-VRC01", "EM4C04-VRC01")
    }
    dftmpstmpVRC01 <- dftmpstmpVRC01[order(match(dftmpstmpVRC01$pair, order_pairs)), ]
    dftmpstmpVRC01$clonotype <- clonotype_info

    dftmpstmpVRC01$pair <- factor(dftmpstmpVRC01$pair, levels=order_pairs)
    dftmpstmpVRC01$id <- c(1:nrow(dftmpstmpVRC01))
    label_data <- dftmpstmpVRC01
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)

    base_data <- dftmpstmpVRC01 %>% 
    group_by(clonotype) %>% 
    summarize(start=min(id), end=max(id) +0.3) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
    base_data <- as.data.frame(base_data)   
    p <- ggplot(dftmpstmpVRC01, aes(x=factor(id), y=avg, fill=avg)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity") +
    ylim(-100,120) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0,1,0,1), "cm") 
    ) + scale_fill_distiller(palette="RdYlGn",breaks=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100), direction = 1) +
    coord_polar(start = 0) + labs(fill="% residual\nbinding") +
    geom_text(data=label_data, aes(x=id, y=avg+10, label=pair, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
    geom_segment(data=base_data, aes(x = start, y = -16, xend = end, yend = -16), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -50, label=clonotype), hjust=c(1,1,0.5,0,0,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
    ggsave(file.path(results_folder, paste0("VRC01_",label,".png")), width=8.5, height=6.5, dpi=300, bg="white")
    ggsave(file.path(results_folder, paste0("VRC01_",label,".pdf")), width=8.5, height=6.5, dpi=300, bg="white")
    
    dftmpstmpVRC01binary <- dftmpstmpVRC01
    dftmpstmpVRC01binary$binary <- dftmpstmpVRC01binary$avg
    binary <- c()
    for (i in dftmpstmpVRC01binary$binary) {
        if ( i > 50) {
            binary <- c(binary, ">50")
        } else { 
            binary <- c(binary, "<50")
        }
    }
    dftmpstmpVRC01binary$binary <- binary
    dftmpstmpVRC01binary$pair <- str_remove_all(dftmpstmpVRC01binary$pair, "-VRC01")
    
    label_data <- dftmpstmpVRC01binary
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)

    base_data <- dftmpstmpVRC01binary %>% 
    group_by(clonotype) %>% 
    summarize(start=min(id), end=max(id) +0.3) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
    base_data <- as.data.frame(base_data)       
    
    p <- ggplot(dftmpstmpVRC01binary, aes(x=factor(id), y=avg, fill=binary)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity") +
    ylim(-100,120) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0,1,0,1), "cm") 
    ) + scale_fill_manual(values=c(">50"="#026464", "<50"="#895129")) +
    coord_polar(start = 0) + labs(fill="% residual\nbinding") +
    geom_text(data=label_data, aes(x=id, y=avg+10, label=pair, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4.5, angle= label_data$angle, inherit.aes = FALSE ) +
    geom_segment(data=base_data, aes(x = start, y = -16, xend = end, yend = -16), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -50, label=clonotype), hjust=c(1,1,0.5,0,0,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
    ggsave(file.path(results_folder, paste0("VRC01_",label,"_binary.png")), width=8.5, height=6.5, dpi=500, bg="white")
    ggsave(file.path(results_folder, paste0("VRC01_",label,"_binary.pdf")), width=8.5, height=6.5, dpi=500, bg="white")

}

count_values<-function(dfcomb, allnormvalues) {
    dfcountvalues <- as.data.frame(matrix(nrow=nrow(dfcomb), ncol=4))
    colnames(dfcountvalues) <- c("combpair", "primary", "secondary", "value.count")
    for (r in rownames(dfcomb)) {
        comb <- dfcomb[r, "combpair"]
        p <- dfcomb[r,"primary"]
        s <- dfcomb[r,"secondary"]
        c = length(allnormvalues[(allnormvalues$primary==p & allnormvalues$secondary==s), "values"])
        dfcountvalues[r, "combpair"] <- comb
        dfcountvalues[r, "primary"] <- p
        dfcountvalues[r, "secondary"] <- s
        dfcountvalues[r, "value.count"] <- c
    }
    write.csv(dfcountvalues, file.path(results_folder, "replicatecount.csv"), quote=F)
    return(dfcountvalues)
}

#############-- Process data --#############

results_folder="1.processedBLIdata_1A8_noRun9-10-11_FULL"
generate_folder(results_folder)

data = read_in_raw_bli_data("./TandemEPComparisons_AM_LWedit_rawdata_FULL.xlsx", results_folder = results_folder)

primary <- unique(data$total_antibodies[!(data$total_antibodies %in% c("SAMPLE"))])
secondary <- unique(data$total_antibodies[!(data$total_antibodies %in% c("SAMPLE"))])
secondary <- secondary[secondary!="EM4C04"]
orderp <- c("70", "4", "12", "18", "1G3", "1A8", "32", "76", "16", "14", "6", "10", "28", "VRC01", "EM4C04")
orders <- c("70", "4", "12", "18", "1G3", "1A8", "32", "76", "16", "14", "6", "10", "28", "VRC01")
dfcomb <- generate_pair_combinations(orderp, orders)

#organize data 
dataorg <- organize_data(data, dfcomb, primary, secondary)
dfcountvalues = count_values(dfcomb, dataorg$allnormvalues)

all_stats <- absolute_diff_between_pairs(dataorg$normavgs, "all")
NoNeg_stats <- absolute_diff_between_pairs(dataorg$normavgsNoNegValues, "NoNegValues")
Zero_stats <- absolute_diff_between_pairs(dataorg$normavgsZeros, "Zeros")

dfall <- generate_initial_avgs_hm(dataorg$normavgs, primary, secondary, "all")
dfNoNeg <- generate_initial_avgs_hm(dataorg$normavgsNoNegValues, primary, secondary, "NoNegValues")
dfZeros <- generate_initial_avgs_hm(dataorg$normavgsZeros, primary, secondary, "Zerovalues")
stop()

generate_final_avgs_hm(dfall, primary, secondary, "all")
generate_final_avgs_hm(dfNoNeg, primary, secondary, "NoNegValues")
generate_final_avgs_hm(dfZeros, primary, secondary, "Zerovalues")

#############-- VRC01 data only --#############

generate_VRC01figs(all_stats, label="all") 
generate_VRC01figs(NoNeg_stats, label="NoNegValues") 
generate_VRC01figs(Zero_stats, label="Zeros") 


