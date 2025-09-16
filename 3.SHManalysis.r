############## --Libraries--##############
options(future.globals.maxSize= 11377049600)
library(ggplot2)
library(Biostrings)
library(seqinr)
library(stringr)
library(dplyr)

############## --Functions--##############
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

theme_Publication <- function(base_size = 14, base_family = "arial") {
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
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}


############## --Run percent identity for clone's of interest SOSIP--##############
sequences_list <- read.fasta(file = "./all_database_Bosinger_LWedit.fa", seqtype = "DNA") # Or "AA" for amino acids

df <- read.csv("./RUp16_results/VDJAnalysisNoWeek36D1_no_42-96/10xinfo&scRepinfo/merged_final_seqs_full.csv", row.names=1)
sosipcells <- read.csv("./RUp16_results/VDJAnalysisNoWeek36D1_no_42-96_final/bar_SOSIP_cells.csv", row.names=1)
df <- df[df$barcode %in% sosipcells$x,]
print(dim(df))

dftotal <- c() 
for (i in rownames(df)) {
    # if (df[i, "sample"]=="RUp16_Wk18") {
        vh <- df[i, "v_gene_heavy"]
        x =str_detect(names(sequences_list), vh)
        y = names(sequences_list)[x]
        yfull = paste(y, collapse=",")
        print(yfull)

        seq2tmp =c(df[i, "fwr1_nt_heavy"], df[i, "cdr1_nt_heavy"], 
            df[i, "fwr2_nt_heavy"], df[i, "cdr2_nt_heavy"], 
            df[i, "fwr3_nt_heavy"])#, df[1, "cdr3_nt_heavy"],  df[1, "fwr4_nt_heavy"]  )
        seq2tmp <- paste(seq2tmp, collapse="")
        pairwiseidentity <- c()
        for (n in y) {
            seqgene <- paste(sequences_list[[n]], collapse="")
            seqgene <- paste(seqgene, collapse="")
            seq1 <- DNAString(seqgene)
            seq2 <- DNAString(seq2tmp)
            alignment <- pwalign::pairwiseAlignment(seq1, seq2, type = "global")
            percent_identity <- pwalign::pid(alignment, type = "PID2")
            pairwiseidentity <- c(pairwiseidentity, percent_identity)
        }

        s <- str_remove_all(df[i, "sample"], "v2")
        s <- str_replace_all(s, "WK", "Wk")
        dftotal <- rbind(dftotal, c(s, max(pairwiseidentity), seqgene, seq2tmp, yfull))
#   }
}
dftotal <- as.data.frame(dftotal)

colnames(dftotal) <- c("sample", "percent_identity", "gene_ref_sequence", "gene_sequence", "gene_names" )
write.csv(dftotal, file.path(results_folder, "datasosip.csv"))
dftotal <- read.csv(file.path(results_folder, "datasosip.csv"))
dftotal$percent_identity <- as.numeric(dftotal$percent_identity)
ggplot(dftotal, aes(x=sample, y=percent_identity)) +
    scale_x_discrete(labels = c("Wk20","Wk26", "Wk36","Wk48","Wk49","Wk80","Wk81","LN-Wk83", "Wk87", "Wk93")) +
    geom_violin() + theme_Publication() + theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9.5),
    axis.title.x=element_blank()) + labs(y="percent identity") +
    stat_summary(fun = mean, geom = "crossbar", width=0.5,color = "black") 

    # stat_summary(fun = median, geom = "crossbar", width=0.5,color = "darkblue")
ggsave(file.path(results_folder, "percentidenitySOSIP.png"), width=5, height=2.5, dpi=300)
ggsave(file.path(results_folder, "percentidenitySOSIP.pdf"), width=5, height=2.5, dpi=300)

df1 <- dftotal %>%
        group_by(sample) %>%
        dplyr::summarize(Mean = mean(percent_identity, na.rm=TRUE),Median = median(percent_identity, na.rm=TRUE))

df1 <- as.data.frame(df1)
write.csv(df1, file.path(results_folder, "datasummarySOSIP.csv"), quote=F,row.names=F)
wildf <- c()
ttestdf <- c()
timepoints <- unique(dftotal$sample)[2:10]
for (t in timepoints) {
    a <- dftotal[dftotal$sample=="RUp16_Wk20", "percent_identity"]
    b  <- dftotal[dftotal$sample==t, "percent_identity"]
   w.stat <-  wilcox.test(a, b)
   t.stat <-  t.test(a, b)
    wildf <- rbind(wildf, c("RUp16_Wk20", t,w.stat$statistic, w.stat$p.value ))
    ttestdf <- rbind(ttestdf, c("RUp16_Wk20", t,t.stat$statistic, t.stat$p.value ))
}
wildf <- as.data.frame(wildf)
ttestdf <- as.data.frame(ttestdf)

colnames(wildf) <- c("timepoint1", "timepoint2", "wilcox.statistic", "p.value")
colnames(ttestdf) <- c("timepoint1", "timepoint2", "t.statistic", "p.value")
wildf$p.val.adj <- p.adjust(wildf$p.value, method="BH")
ttestdf$p.val.adj <- p.adjust(ttestdf$p.value, method="BH")
write.csv(wildf, file.path(results_folder, "wilcoxstatisticsSOSIP.csv"), quote=F)
write.csv(ttestdf, file.path(results_folder, "tstatisticsSOSIP.csv"), quote=F)

############## --Run percent identity for clone's neutralizing clonotypes of interest  --##############
sequences_list <- read.fasta(file = "./all_database_Bosinger_LWedit.fa", seqtype = "DNA") # Or "AA" for amino acids

df <- read.csv("./RUp16_results/VDJAnalysisNoWeek36D1_no_42-96/10xinfo&scRepinfo/merged_final_seqs_full.csv", row.names=1)
coi <- c("684", "695", "1303","386","385","599")
df <- df[df$LW_clones %in% coi,]
print(dim(df))
dftotal <- c() 
for (i in rownames(df)) {
    # if (df[i, "sample"]=="RUp16_Wk18") {
        vh <- df[i, "v_gene_heavy"]
        x =str_detect(names(sequences_list), vh)
        y = names(sequences_list)[x]
        yfull = paste(y, collapse=",")
        print(yfull)

        seq2tmp =c(df[i, "fwr1_nt_heavy"], df[i, "cdr1_nt_heavy"], 
            df[i, "fwr2_nt_heavy"], df[i, "cdr2_nt_heavy"], 
            df[i, "fwr3_nt_heavy"])#, df[1, "cdr3_nt_heavy"],  df[1, "fwr4_nt_heavy"]  )
        seq2tmp <- paste(seq2tmp, collapse="")
        pairwiseidentity <- c()
        for (n in y) {
            seqgene <- paste(sequences_list[[n]], collapse="")
            seqgene <- paste(seqgene, collapse="")
            seq1 <- DNAString(seqgene)
            seq2 <- DNAString(seq2tmp)
            alignment <- pwalign::pairwiseAlignment(seq1, seq2, type = "global")
            percent_identity <- pwalign::pid(alignment, type = "PID2")
            pairwiseidentity <- c(pairwiseidentity, percent_identity)
        }

        s <- str_remove_all(df[i, "sample"], "v2")
        s <- str_replace_all(s, "WK", "Wk")
        dftotal <- rbind(dftotal, c(s, max(pairwiseidentity), seqgene, seq2tmp, yfull, df[i, "LW_clones"]))
#   }
}
dftotal <- as.data.frame(dftotal)

colnames(dftotal) <- c("sample", "percent_identity", "gene_ref_sequence", "gene_sequence", "gene_names", "clonotype" )
write.csv(dftotal, file.path(results_folder, "datacoi.csv"))
dftotal <- read.csv(file.path(results_folder, "datacoi.csv"))
dftotal$percent_identity <- as.numeric(dftotal$percent_identity)
coi <- c("684","386", "695","385","1303","599")

dftotal$clonotype <- factor(dftotal$clonotype, levels=coi)
ggplot(dftotal, aes(x=sample, y=percent_identity)) +
    scale_x_discrete(labels = c("Wk20", "Wk26", "Wk36","Wk48","Wk49","Wk80","Wk81","LN-Wk83", "Wk87", "Wk93")) +
    geom_violin() + theme_Publication() + theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9.5),
    axis.title.x=element_blank()) + labs(y="percent identity") + facet_wrap(~clonotype, ncol=2) +
    stat_summary(fun = mean, geom = "crossbar", width=0.5,color = "black") 
ggsave(file.path(results_folder, "percentidenitycoi.png"), width=8, height=6.5, dpi=300)
ggsave(file.path(results_folder, "percentidenitycoi.pdf"), width=8, height=6.5, dpi=300)

df1 <- dftotal %>%
        group_by(clonotype,sample ) %>%
        dplyr::summarize(Mean = mean(percent_identity, na.rm=TRUE),Median = median(percent_identity, na.rm=TRUE))

df1 <- as.data.frame(df1)
write.csv(df1, file.path(results_folder, "datasummarycoi.csv"), quote=F,row.names=F)
wildf <- c()
ttestdf <- c()
for (c in unique(as.character(dftotal$clonotype))) {
    dftotalsub <- dftotal[dftotal$clonotype==c,]
    if (c %in% c("684", "695")) {
        timepoints <- unique(dftotal$sample)[3:10]
        for (t in timepoints) {
            a <- dftotalsub[dftotalsub$sample=="RUp16_Wk26", "percent_identity"]
            b  <- dftotalsub[dftotalsub$sample==t, "percent_identity"]
            w.stat <-  wilcox.test(a, b)
            # t.stat <-  t.test(a, b)
            wildf <- rbind(wildf, c("RUp16_Wk26", t,w.stat$statistic, w.stat$p.value, c ))
            # ttestdf <- rbind(ttestdf, c("RUp16_Wk26", t,t.stat$statistic, t.stat$p.value, c ))
        }
    } else if (c %in% c("386", "385")) {
        timepoints <- unique(dftotal$sample)[2:10]
        for (t in timepoints) {
            a <- dftotalsub[dftotalsub$sample=="RUp16_Wk20", "percent_identity"]
            b  <- dftotalsub[dftotalsub$sample==t, "percent_identity"]
            w.stat <-  wilcox.test(a, b)
            # t.stat <-  t.test(a, b)
            wildf <- rbind(wildf, c("RUp16_Wk20", t,w.stat$statistic, w.stat$p.value, c ))
            # ttestdf <- rbind(ttestdf, c("RUp16_Wk26", t,t.stat$statistic, t.stat$p.value, c ))
        }
    } else if (c %in% c("599")) {
        timepoints <- unique(dftotal$sample)[c(3:5, 7:8,10)]
        for (t in timepoints) {
            a <- dftotalsub[dftotalsub$sample=="RUp16_Wk26", "percent_identity"]
            b  <- dftotalsub[dftotalsub$sample==t, "percent_identity"]
            w.stat <-  wilcox.test(a, b)
            # t.stat <-  t.test(a, b)
            wildf <- rbind(wildf, c("RUp16_Wk20", t,w.stat$statistic, w.stat$p.value, c ))
            # ttestdf <- rbind(ttestdf, c("RUp16_Wk26", t,t.stat$statistic, t.stat$p.value, c ))
        }
    } else {
        timepoints <- unique(dftotal$sample)[c(5, 7:10)]
        for (t in timepoints) {
            a <- dftotalsub[dftotalsub$sample=="RUp16_Wk48", "percent_identity"]
            b  <- dftotalsub[dftotalsub$sample==t, "percent_identity"]
            w.stat <-  wilcox.test(a, b)
            # t.stat <-  t.test(a, b)
            wildf <- rbind(wildf, c("RUp16_Wk48", t,w.stat$statistic, w.stat$p.value, c ))
            # ttestdf <- rbind(ttestdf, c("RUp16_Wk8", t,t.stat$statistic, t.stat$p.value, c ))
        }
    }
   
}
wildf <- as.data.frame(wildf)
# ttestdf <- as.data.frame(ttestdf)

colnames(wildf) <- c("timepoint1", "timepoint2", "wilcox.statistic", "p.value", "clonotype")
# colnames(ttestdf) <- c("timepoint1", "timepoint2", "t.statistic", "p.value", "clonotype")
wildf$p.val.adj <- p.adjust(wildf$p.value, method="BH")
# ttestdf$p.val.adj <- p.adjust(ttestdf$p.value, method="BH")
write.csv(wildf, file.path(results_folder, "wilcoxstatisticscoi.csv"), quote=F)
# write.csv(ttestdf, file.path(results_folder, "tstatisticscoi.csv"), quote=F)

neut <- c()
for (c in dftotal$clonotype) {
    if (c %in% c("684","695", "1303")) {
        neut <- c(neut, "Neutralizer")
    } else {
        neut <- c(neut, "NonNeutralizer")
    }
}
dftotal$neut <- neut
ggplot(dftotal, aes(x=sample, y=percent_identity, color=neut)) +
    scale_x_discrete(labels = c("Wk20", "Wk26", "Wk36","Wk48","Wk49","Wk80","Wk81","LN-Wk83", "Wk87", "Wk93")) +
    scale_color_manual(values=c("Neutralizer"="#e09100", "NonNeutralizer"="#690da3")) +
    geom_violin(position=position_dodge(1)) + theme_Publication() + theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9.5),
    axis.title.x=element_blank(), legend.position = "top", legend.direction = "horizontal", legend.title=element_blank()) + labs(y="percent identity") + 
    stat_summary(fun = mean, geom = "crossbar", width=0.5,position=position_dodge(1),show.legend = FALSE) 
ggsave(file.path(results_folder, "percentidenitycoineut.png"), width=5.5, height=3, dpi=300)
ggsave(file.path(results_folder, "percentidenitycoineut.pdf"), width=5.5, height=3, dpi=300)
wildf <-c()
timepoints <- unique(dftotal$sample)[2:10]
for (t in timepoints) {
    a <- dftotal[(dftotal$sample==t & dftotal$neut=="Neutralizer"), "percent_identity"]
    b  <- dftotal[(dftotal$sample==t & dftotal$neut=="NonNeutralizer"), "percent_identity"]
    w.stat <-  wilcox.test(a, b)
    wildf <- rbind(wildf, c(t,w.stat$statistic, w.stat$p.value))
}

wildf <- as.data.frame(wildf)
colnames(wildf) <- c("timepoint", "wilcox.statistic", "p.value")
wildf$p.val.adj <- p.adjust(wildf$p.value, method="BH")
write.csv(wildf, file.path(results_folder, "wilcoxstatisticNeut.csv"), quote=F)