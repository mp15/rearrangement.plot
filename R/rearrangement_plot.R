# An R file to be sourced.
#
# Rearrangement plot (aka chromothripsis paper type plot) for plotting 
# rearrangement breakpoints and copy number.
#
# After sourcing workspace will have a variable 'chr_lens', and three functions:
# arc(), window_means() and plot_rearrangements()
# 
# Author: Yilong Li
###############################################################################

# To read in this file: 
# source("/Users/yl3/Documents/workspace/rg_ordering_pilot/rearrangement_plot.R")

# /// EDITED by Simon Brunner /// # 


library(dplyr)
library(scales)
library(quantsmooth)  # For ideogram

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=scales::alpha)) 
}

arc = function(x0, x1, y, xr, yr, col, lwd) {
    x = (x0 + x1)/2  # Center of arc
    xr = x - x0 	 # x-radius of arc
    
    apply(
        cbind(x, y, xr, yr, col),
        1,
        function(z) {
            x   = as.numeric(z[1])
            y   = as.numeric(z[2])
            xr  = as.numeric(z[3])
            yr  = as.numeric(z[4])
            col = z[5]
            x_points = seq(x - xr, x + xr, length.out = 200)
            y_points = y + yr * sqrt( 1  -  ( (x_points-x) / xr )^2 )
            
            lines(
                x_points,
                y_points,
                col = col,
                lwd = lwd
            )
        }
    )
    
    return()
}

window_means = function (coords, cns, min_pos, max_pos, win_size) {
    cut_levels = cut(
        coords,
        seq(min_pos, max_pos + win_size, win_size),
        labels = F,
        include_lowest = T,
        right = F
    )
    
    cut_values = by(
        cns,
        cut_levels,
        mean
    )
    
    return(
        cut_values[ as.character(1:ceiling((max_pos-min_pos+1)/win_size)) ]
    )
}

#' Plot rearrangements from SV and CNV data
#'
#' This function plots structural variant and CNV data against the genome..
#' @param bedpe File handle for BEDPE containing SVs
#' @param chrs The chromosomes to plot for
#' @param chr_lens File handle for file containing chromosome_lengths
#' @param cn_bedgraph File handle for bg file contain copy number information. Defaults to NULL. 
#' @param CN_SV_gap: (default: False) If set to true then a gap of size yrange_size will be left between copy number (CN) and structural variant (SV) data.
#' @param specific_SV_cols (default: NA): if set, then SVs will be plotted with the colors provided in the vector (one entry for each SV)
#' @param y_sv1: yrange[2]*(1+y_sv1) is where the lower SVs are plotted
#' @param y_sv2: yrange[2]*(1+y_sv2) is where the higher SVs are plotted
#' @param yaxis_side: side where to plot y axis. Takes same values as side argument in mtext (1=bottom, 2=left, 3=top, 4=right)
#' @keywords rearrangment
#' @export
#' @examples
#' fn_bedpe <- '/data/CGP_donor_A/sanger/A_SV_4_0_12/A_vs_B.annot.bedpe'
#' fn_cn <- '/data/CGP_donor_A/sanger/A_SV_4_0_12/intermediates/B.ngscn.abs_cn.bg.gz'
#' bedpe <- read.table(fn_bedpe, header = F, sep = "\t", stringsAsFactors = F)
#' zz <- gzfile('file.csv.gz','rt')  
#' cn <- read.table(zz,    header = F, sep = "\t", stringsAsFactors = F, colClasses = c("factor", "numeric", "numeric", "numeric"))
#' chr_lens <- read.table('~/dev/params/hg19.chrom_sizes.txt', header=F, sep="\t", row.names=1, colClasses = c("character","character", "numeric"))
#' temp <- rownames(chr_lens)
#' chr_lens <- chr_lens[,2]
#' names(chr_lens) <- temp
#' this_sample <- 'test'
#' fpath_out <- '~/dev/testplot'
#' pdf(paste0(fpath_out,".pdf"), h = 4, w = 10)
#' plot_rearrangements(bedpe = bedpe, chrs=c(1), chr_lens = chr_lens, cn_bedgraph = cn, cn_win_size = 1e4, cn_cex = 0.3)
#' dev.off()

plot_rearrangements = function(
    bedpe, chrs, chr_lens, cn_bedgraph = NULL, segments = NULL,
    yrange = NULL, ideogram=T, cn_cex=0.5, lwd = 0.75, cn_win_size = 1e5,
    BFB_ids = c(), arrow_ln = 0.15, xlim = NULL, chr_lim = NULL, annot = NULL, main = NULL, CN_SV_gap=F,
    ascat_tbl=NULL, specific_SV_cols=NA, p_ylim_quantile=0.999, yaxis_ticks=NA, y_sv1 = 0.3, y_sv2 = 0.75,
    yaxis_side = 4
) {
  # Args:
    chrs = as.character(chrs)
    if (!is.null(chr_lim)) {
        chr_lim = as.character(chr_lim)
    }
    chr_cum_lns = c(0, cumsum(chr_lens[chrs])[-length(chrs)])
    names(chr_cum_lns) = chrs
    xrange_size = cumsum(chr_lens[chrs])
    
    if (is.null(xlim)) {
        xlim = c(1, sum(chr_lens[chrs]))
        
        if (is.null(yrange)) {
            yrange = c(0, quantile(cn_bedgraph[cn_bedgraph[,1] %in% chrs, 4], p = p_ylim_quantile))
            yrange = c(floor(yrange[1]), ceiling(yrange[2]))
            if(yrange[2]<2) {
              yrange[2] = 2
            }
        }
    } else if (length(chrs) > 1) {
        if (!is.null(chr_lim)) {
            stop()
        }
        
        if (is.null(yrange)) {
            yrange = c(0, quantile(cn_bedgraph[cn_bedgraph[,1] %in% chr_lim & cn_bedgraph[,2] > xlim[1] & cn_bedgraph[,3] < xlim[2], 4], p = p_ylim_quantile))
            yrange = c(floor(yrange[1]), ceiling(yrange[2]))
        }
        
        xlim = chr_cum_lns[chr_lim] + xlim
    } else {
        if (is.null(yrange)) {
            yrange = c(0, quantile(cn_bedgraph[cn_bedgraph[,1] == chrs & cn_bedgraph[,2] > xlim[1] & cn_bedgraph[,3] < xlim[2], 4], p = p_ylim_quantile))
            yrange = c(floor(yrange[1]), ceiling(yrange[2]))
        }
    }
    
    # If CN_SV_gap is False, then the y-axis positioning of SV calls refers to the default yrange
    if(CN_SV_gap==F) {
      yrange_CN = yrange
      yrange_CN[1] = yrange_CN[1]-0.5
    # If CN_SV_gap is True, then the y-axis positioning of SV calls refers to a different yrange, which extends twice the size of the standard yrange
    } else if(CN_SV_gap==T) {
      yrange_CN = yrange
      yrange_CN[2] = yrange_CN[2] + yrange[2]
    }
    
    if (!is.null(cn_bedgraph)) {
        cn = cn_bedgraph[cn_bedgraph[,1] %in% chrs, ]
    } else if (!is.null(segments)) {
    } else {
        stop("Either cn_bedgraph or segments must be provided")
    }
    
    yrange_size = yrange[2] - yrange[1]
    yrange_size_CN = yrange_CN[2] - yrange_CN[1]
    
    td_col         = "darkorange4"
    del_col        = "darkslateblue"
    inter_chrs_col = "darkorchid1"
    tail_tail_col  = "cadetblue4"
    head_head_col  = "chartreuse2"
    
    # Determine the SV class based on BRASS2 strand calls, where:
    # D: deletion; TD: tandem-duplication; HH: head-to-head inversion; TT: tail-to-tail inversion
    bedpe = dplyr::mutate(bedpe, SV_class=ifelse(V12=='translocation', 'translocation',
                                   ifelse(V9=='+' & V10=='+', 'D', 
                                   ifelse(V9=='-' & V10=='-', 'TD', 
                                   ifelse(V9=='+' & V10=='-', 'HH','TT')))))
    
    # Create the plot
    # par(mar = c(5, 4, 2, 2) + .1)
    
    # Define ylim depending on whether ascat alleles are provided
    if(is.null(ascat_tbl)) {
      def_ylim = c(yrange[1] - .2*yrange_size, yrange_CN[2] + 1.4*yrange_size)
    } else {
      def_ylim = c(yrange[1] - .4*yrange_size, yrange_CN[2] + 1.4*yrange_size)
    }
    
    plot(
        c(),
        ylim = def_ylim,
        xlim = xlim,
        bty  = "n",
        yaxt = "n",
        xaxt = "n",
        xlab = "",
        ylab = "",
        yaxs = "i",
        xaxs = "i",
        main = main
    )
    
    # X axis names and ticks
    par(mgp = par("mgp") + c(0,1,0))
    if (all(xlim == c(1, sum(chr_lens[chrs])))) {
        axis(
            1,
            at = (cumsum(chr_lens[chrs]) + chr_cum_lns)/2,
            labels = paste("chr", chrs, sep=""),
            tick = F,
            cex.lab = 1.5,
            cex.axis = 0.8
        )
    } else if (!is.null(chr_lim)) {
        axis(
            1,
            at = mean(xlim),
            labels = paste("chr", chr_lim, " position (Mb) 2", sep=""),
            tick = F,
            cex.lab = 1.5
        )
    } else {
        axis(
            1,
            at = mean(xlim),
            labels = paste("chr", chrs, " position (Mb) 3", sep=""),
            tick = F,
            cex.lab = 1.5
        )
    }
    par(mgp = par("mgp") - c(0,1,0))
      
    if(CN_SV_gap==T) {
      #par(oma = c(0,0,0,2))
    }
    
    if (length(chrs) > 1) {
        if (all(xlim == c(1, sum(chr_lens[chrs])))) {
            for (c in chrs) {
                pretty_ticks = pretty(c(1, chr_lens[c]))
                pretty_ticks = pretty_ticks[which(pretty_ticks < chr_lens[c])]
                axis(1, at = pretty_ticks + chr_cum_lns[c], labels = pretty_ticks/1e6)
            }
        } else {
            if (is.null(chr_lim) || !(chr_lim %in% names(chr_lens))) {
                stop()
            }
            
            pretty_ticks = pretty(xlim - chr_cum_lns[chr_lim])
            
            axis(1, at = pretty_ticks + chr_cum_lns[chr_lim], labels = pretty_ticks/1e6)
        }
    } else {
        if (all(xlim == c(1, sum(chr_lens[chrs])))) {
            axis(1, at = axisTicks(usr=c(1, chr_lens[chrs]), log=F), labels = axisTicks(usr=c(1, chr_lens[chrs]), log=F)/1e6)
        } else {
            pretty_ticks = pretty(xlim)
            axis(1, at = pretty_ticks, labels = pretty_ticks / 1e6)
        }
    }
    
    
    # Shaded grid
    for (i in yrange[1]:yrange[2]) {
        polygon(
            c(1, sum(chr_lens[chrs]), sum(chr_lens[chrs]), 1),
            i - 0.5 + c(0, 0, 1, 1),
            col=rgb(.1, .1, .1, ifelse(i %% 2 == 0, 0.1, 0.05)),
            lty=0
        )
    }
    
    
    # Line to separate chromosomes
    if (length(chrs) > 1) {
        segments(
            x0 = cumsum(chr_lens[chrs])[-length(chrs)],
            y0 = yrange[1] - 0.5,
            y1 = yrange[2] + 0.5
        )
    }
    
    # If provided, plot ASCAT allele data
    if(!is.null(ascat_tbl)) {
      #points(x=ascat_tbl$Position, y=(0.2*yrange_size*ascat_tbl$BAF)-0.25*yrange_size, col=rgb(0,0,0,0.05), pch=16, cex=0.2)
      points(x=ascat_tbl$Position, y=(0.2*yrange_size*ascat_tbl$BAF)-0.25*yrange_size, col=add.alpha(ascat_tbl$col, alpha=0.5), pch=16, cex=0.2)
      #axis(side = 4, at=c(-0.25*yrange_size, -0.05*yrange_size), labels=c(0,1), las=2, lwd=0.5, cex.axis=1.0)
      mtext(side=2, 'Allele',at=-0.125*yrange_size, cex=1.0, las=2, line=0.2)
    }
    
    yaxis_line = ifelse(yaxis_side == 2, 2.2, 0.2)
    if(CN_SV_gap) {
      print(par('cex.axis'))
      #mtext(side=2, 'CN',at=yrange_size/2, line=3, cex=par('cex.axis'), las=2)
      mtext(side=2, 'Copy\nnum.',at=yrange_size/2, cex=1.0, las=2, line=yaxis_line)
    } else {
      #title(ylab = "CN", line=3, cex=par('cex.axis'), las=2)
      mtext(side=2, 'Copy\nnum.',at=yrange_size/2, cex=1.0, las=2, line=yaxis_line)
    }
    
    # Plot rearrangements: First dotted lines
    segments(
        x0 = c(1, 1),
        x1 = c(sum(chr_lens[chrs]), sum(chr_lens[chrs])),
        y0 = c(yrange_CN[2] + y_sv1*yrange_size, yrange_CN[2] + y_sv2*yrange_size),
        lty = 3
    )
    abline(h = yrange[2] + 0.5)
    
    # Then 'intra-chromosomal' translocations
    sel = bedpe[,1] %in% chrs & bedpe[,4] %in% chrs & !(bedpe[,7] %in% BFB_ids)  # Breakage-fusion-bridge RGMTS to be plotted separately
    if(is.na(specific_SV_cols)) {
      col =
        ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
                ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
                        ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col, tail_tail_col)))
    } else {
      col = specific_SV_cols[sel]
    }
    
    if (sum(sel) > 0) {
      arc(
        x0  = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
        x1  = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
        y   = ifelse(bedpe[sel,9] == bedpe[sel, 10], yrange[2]+.3*yrange_size, yrange[2]+.75*yrange_size),
        yr  = ifelse(bedpe[sel, 10] == "-", 1, -1) * 0.2 * yrange_size,
        col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col},
        lwd = lwd
      )
      segments(
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
        y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75*yrange_size, 0.3*yrange_size),
        y1 = yrange[1],
        col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col},
        lwd = lwd
      )
      segments(
        x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
        y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75*yrange_size, 0.3*yrange_size),
        y1 = yrange[1],
        col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col},
        lwd = lwd
      )
    }
    #----
    
    # Then rearrangments where low end %in% chrs and !(high end %in% chrs) 
    sel = bedpe[,1] %in% chrs & !(bedpe[,4] %in% chrs)
    if (sum(sel) > 0) {
        # arrows(
    	if(is.na(specific_SV_cols)) {
            col = rgb(t(col2rgb("black")), alpha = 127, max=255)
		} else {
			col = specific_SV_cols[sel]
		} 
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange_CN[1],
            y1 = yrange_CN[2] + yrange_size,
            col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col},
            #col = rgb(t(col2rgb("black")), alpha = 127, max=255),
            lwd = lwd  # ,
            # length = arrow_ln
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange_CN[2] + yrange_size,
            x1 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y1 = ifelse(bedpe[sel, 9] == "+", yrange_CN[2] + 1.1 * yrange_size, yrange_CN[2] + 1.3 * yrange_size),
            xpd = NA
        )
        text(
            as.character(bedpe[sel, 4]),
            x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            # y = yrange[2] + 1.2 * yrange_size
            y = yrange_CN[2] + ifelse(bedpe[sel, 9] == "+", 1.2, 1.4) * yrange_size,
            xpd = NA, cex = par('cex.axis')
        )
    }
    
    # Then rearrangments where high end %in% chrs and !(low end %in% chrs) 
    sel = !(bedpe[,1] %in% chrs) & bedpe[,4] %in% chrs
    if (sum(sel) > 0) {
        # arrows(
    	if(is.na(specific_SV_cols)) {
            col = rgb(t(col2rgb("black")), alpha = 127, max=255)
		} else {
			col = specific_SV_cols[sel]
		} 
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y0 = yrange_CN[1],
            y1 = yrange_CN[2] + yrange_size,
            col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col},
            lwd = lwd  # ,
            # length = arrow_ln
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y0 = yrange_CN[2] + yrange_size,
            x1 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y1 = ifelse(bedpe[sel, 10] == "+", yrange_CN[2] + 1.1 * yrange_size, yrange_CN[2] + 1.3 * yrange_size),
            xpd = NA
        )
        text(
            as.character(bedpe[sel, 1]),
            x = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            # y = yrange[2] + 1.2 * yrange_size
            y = yrange_CN[2] + ifelse(bedpe[sel, 10] == "+", 1.2, 1.4) * yrange_size,
            xpd = NA, cex = par('cex.axis')
        )
    }
    
    
    # Then BFBs, currently only supporting 2 BFB events max
    # Also, in case of BFBs, only plotting 1st end
    if (length(BFB_ids) > 2) {
        stop()
    }
    if (length(BFB_ids) == 2) {
        sel = bedpe[,7] == BFB_ids[2]
    	if(is.na(specific_SV_cols)) {
            col = ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
                ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
                    ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col, 		 tail_tail_col)))
		} else {
			col = specific_SV_cols[sel]
		} 
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange_CN[2] + 1.2*yrange_size,
            y1 = yrange_CN[1],
            col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col}
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange_CN[2] + 1.2*yrange_size,
            y1 = yrange_CN[2] + (1.2 + 0.1)*yrange_size,
            col = if(is.na(specific_SV_cols)) {rgb(t(col2rgb(col)), alpha = 127, max=255)} else {col}
        )
        
        # The curved arrow
        lines(
            x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(seq(-pi/2, pi/2, length.out=20)) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50,
            y = yrange_CN[2] + (1 + 0.3)*yrange_size +  0.05*yrange_size*sin(seq(-pi/2, pi/2, length.out=20)),
            col = col,
            lwd = 2 * lwd
        )
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50
        arrows(
            x0 = x0,
            x1 = x0 + ifelse(bedpe[sel,9] == "+", -1, 1) * xrange_size / 50,
            y0 = yrange_CN[2] + (1 + 0.3)*yrange_size +  0.05*yrange_size*sin(-pi/2),
            col = col,
            lwd = 2 * lwd,
            length = arrow_ln
        )
    }
    if (length(BFB_ids) >= 1) {
        abline(h = yrange_CN[2] + c(1, 1.2)*yrange_size)
        if (length(BFB_ids) > 1) abline(h = yrange_CN[2] + 1.4*yrange_size)
        
        sel = bedpe[,7] == BFB_ids[1]
        col = 
            ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
                ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
                    ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col, 		 tail_tail_col)))
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange_CN[2] + 1*yrange_size,
            y1 = yrange_CN[1],
            col = rgb(t(col2rgb(col)), alpha = 127, max=255)
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange_CN[2] + 1*yrange_size,
            y1 = yrange_CN[2] + (1 + 0.1)*yrange_size,
            col = col
        )
        
        # The curved arrow
        lines(
            x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(seq(-pi/2, pi/2, length.out=20)) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50,
            y = yrange_CN[2] + (1 + 0.1)*yrange_size +  0.05*yrange_size*sin(seq(-pi/2, pi/2, length.out=20)),
            col = col,
            lwd = 2 * lwd
        )
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50
        arrows(
            x0 = x0,
            x1 = x0 + ifelse(bedpe[sel,9] == "+", -1, 1) * xrange_size / 50,
            y0 = yrange_CN[2] + (1 + 0.1)*yrange_size +  0.05*yrange_size*sin(-pi/2),
            col = col,
            lwd = 2 * lwd,
            length = arrow_ln
        )
    }
    
    
    # Annotations on top
    if (!is.null(annot) && sum(annot[,1] %in% chrs > 0)) {
        sel = annot[,1] %in% chrs
        segments(
            x0 = chr_cum_lns[as.character(annot[sel, 1])] + annot[sel, 2],
            y0 = 1.1 * yrange_size + yrange_CN[2] + (0:(sum(annot[,1] %in% chrs)-1)) * 0.05 * yrange_size,
            x1 = chr_cum_lns[as.character(annot[sel, 1])] + annot[sel, 3],
            lwd = 2
        )
        
        text(
            chr_cum_lns[as.character(annot[sel, 1])] + rowMeans(annot[sel, 2:3]),
            1.2 * yrange_size + yrange_CN[2] + (0:(sum(annot[,1] %in% chrs)-1)) * 0.05 * yrange_size,
            labels = annot[sel,4], cex = par('cex.axis')
        )
    }

    # Finally ideogram
    if (xlim[2] - xlim[1] < 10e6) {
        print("Ideogram plotting disabled because xlim[2] - xlim[1] < 10e6")
        
        ideogram = F
    }
    if (ideogram) {
        for (c in chrs) {
            quantsmooth::paintCytobands(
                c,
                pos = c(1 + chr_cum_lns[c], yrange[1]-ifelse(is.null(ascat_tbl),0.1,0.3)*yrange_size),
                units="bases",
                width = 0.1*yrange_size,
                length.out = chr_lens[c],
                legend = F
            )
        }
    }
    
    # Plot CN
    win_size = cn_win_size
    for (c in chrs) {
      if (!is.null(cn_bedgraph)) {
        # sel = cn[,1] == c
        # When plotting only one chromome with zoomed-in image, only plot
        # the data that will be visible in the graph.
        if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
          sel = cn[,1] == c & cn[,2] > xlim[1] - win_size & cn[,3] < xlim[2] + win_size
          x = win_size/2 + seq(xlim[1]-win_size, xlim[2]+win_size, win_size)
          y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], xlim[1]-win_size, xlim[2]+win_size, win_size)
          points(
            x = x,
            y = y,
            pch = 16,
            cex = cn_cex,
            xpd=T,
            col=scales::alpha('black', 1)
          )
        }
        else {
          sel = cn[,1] == c
          points(
            x = chr_cum_lns[c] + win_size/2 + seq(1, chr_lens[c], win_size),
            y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], 1, chr_lens[c], win_size),
            pch = 16,
            cex = cn_cex,
            col=scales::alpha('black', 1)
          )
        }
      }
      
      
      
      if (!is.null(segments)) {
        sel = segments[,1] == c
        
        segments(
          x0 = chr_cum_lns[c] + segments[sel, 2] + 1,
          x1 = chr_cum_lns[c] + segments[sel, 3],
          y0 = segments[sel, 4],
          lwd = 2,
          col = "blue"
        )
      }
    }
    # Define Y axis, if exact ticks have been provided
    if(!is.na(yaxis_ticks)) {
      axis(
        4,
        at = yaxis_ticks,
        tick = T,
        las=2
      )
    } else {
      #axis(2, at = axisTicks(usr=yrange, log=F, nint=4), las=2)
      #axis(2, at = pretty(x=c(yrange[1], 0.75*yrange[2]), n=3), las=2, cex.axis=1.3)
      yaxis_ticks = round(axisTicks(usr=c(yrange[1], 0.9*yrange[2]), nint=3, log=F))
      axis(yaxis_side, at = yaxis_ticks, las=2, cex.axis=1)
    }
    
    c(yrange=yrange,yrange_size=yrange_size)
}

