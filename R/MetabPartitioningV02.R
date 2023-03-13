## Neha experiments. Partitioning calculations
# 24 hr exposure time

# Install packages
install.packages("readxl")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("extrafont")

# load libraries
library(readxl)
library(reshape2)
library(ggplot2)
library(extrafont)

# Read data ---------------------------------------------------------------

# Read data.xlsx
d <- data.frame(read_excel("data.xlsx", sheet = "dataV2",
                           col_names = TRUE, col_types = NULL))

# Calculations ------------------------------------------------------------
# Name parameters
congener <- d$congener
logKa.w <- d$logKair.water
dUaw <- d$dUaw
logKlip.w <- d$logKlipid.water # storage lipid
logKpro.w <- d$logKprotein.water # muscle protein

#Fixed parameters
# 24-well plate https://www.corning.com/catalog/cls/documents/drawings/MD_Microplate_Dimension_Sheets_Multiple_Well.pdf
# Multiwell dimensions
Vt <- 3.47/1000 # Total volume L
Vm <- 0.5/1000 # Medium volume L
Va <- Vt-Vm # Air volume L
A <- 1.93 # Area cm2
ht <- 1.74 # Height cm
r <- sqrt((A/pi)) # radius cm
hm <- Vm*1000/(pi*r^2) # Medium height cm

# Function to calculate fractions
fraction = function(logKa.w, dUaw, logKlip.w, logKpro.w) {
  
  # Concentrations in the wells
  Cell <- 4*10^4 # # cell/well
  C.cell <- Cell/Vm # Concentration of cell per well #cell/Vm (L)
  C.lipid.cell <- 9.57*10^-5/10^9 # lipid content kg/#cell
  C.lipid.L <- C.lipid.cell*C.cell # kg/L
  dlipids <- 0.905 # kg/L ask! ref: https://pubmed.ncbi.nlm.nih.gov/8148928/
  C.lipid <- C.lipid.L/dlipids # Lipid concentration Llipid/lwater
  C.prot.cell.0 <- 1.18*10^-4/10^9 # prot content kg/#cell
  C.prot.L <- C.prot.cell.0*C.cell # kg/L
  dprot <- 1.43 # kg/L ask! ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot <- C.prot.L/dprot # Lipid concentration Lprot/lwater
  C.water.cell <- 2.84*10^-6 # uL/cell
  V.water.cell <- C.water.cell*Cell/10^6 # L water inside cell per well
  
  # Temperature correction for Kaw
  R <- 8.3144 # J/(mol K)
  Tst <- 25 # Standard temperature C
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  T2 <- 37 # Incubation temperature C
  T2 <- 273.15 + T2
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/T2-1/Tst.1)) # Ka.w corrected by water and air temps
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculation
  #C.t <- C.dis*(Vm + V.water.cell + 10^(log.Klip.w)*C.lipid + 10^(log.Kpro.w)*C.prot+
  #                10^(log.Kw.air)*Va/(Vm)/(Vm + V.water.cell)
  den <- 2 + 10^(logKlip.w)*C.lipid + 10^(logKpro.w)*C.prot+
    10^(logKa.w.t)*Va/(Vm)
  f.dis <- 1/den
  f.dis.m <- f.dis*(Vm - V.water.cell)/Vm
  f.dis.c <- f.dis*V.water.cell/(Vm)
  f.lip <- 10^(logKlip.w)*C.lipid/den
  f.prot <- 10^(logKpro.w)*C.prot/den
  f.cell <- f.dis.c + f.lip + f.prot
  f.prot.cell <- f.prot/f.cell
  f.lip.cell <- f.lip/f.cell
  f.dis.cell <- f.dis.c/f.cell
  f.air <- 10^(logKa.w.t)*Va/(Vm)/den
  f.well <- f.cell + f.dis.m + f.air
  f.cell.well <- f.cell/f.well
  f.air.well <- f.air/f.well
  f.dis.well <- f.dis.m/f.well
  frac <- c(f.dis.well, f.cell.well, f.air.well, f.prot.cell,
            f.lip.cell, f.dis.cell)
}

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener) {
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKlip.w[i], logKpro.w[i]))
}

final.result <- data.frame(congener, result)
names(final.result) <- c("congener", "fract.dis.well", "fract.cell.well",
                         "fract.air.well", "fract.prot.cell", "fract.lip.cell",
                         "fract.dis.cell")

# Plots -------------------------------------------------------------------
# Transform data.frame p.1 to 3 column data.frame
plot.data <- melt(final.result, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")
# Select only well fractions
plot.data.well <- plot.data[plot.data$phase %in% c("fract.air.well",
                                                   "fract.dis.well",
                                                   "fract.cell.well"), ]

# Organize compounds name
plot.data.well$congener <- factor(plot.data.well$congener,
                            levels = c('PCB3', '4-OH-PCB3', '4-OH-PCB3 sulfate',
                                       'PCB11', '4-OH-PCB11', '4-OH-PCB11 sulfate',
                                       'PCB25', '4-OH-PCB25', '4-OH-PCB25 sulfate',
                                       'PCB52', '4-OH-PCB52', '4-OH-PCB52 sulfate'))
# Organize phases
plot.data.well$phase <- factor(plot.data.well$phase,
                        levels = c('fract.air.well', 'fract.dis.well',
                                   'fract.cell.well'))

ggplot(plot.data.well, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("air" , "medium", "cell"),
                    values = c("deepskyblue", "brown1", "coral4")) +
  theme_classic() +
  theme(aspect.ratio = 10/10,
        text = element_text(size = 14, family = "Helvetica",
                            face = "bold", color = "black")) +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# (2) Just fractions inside cell, i.e., protein, lipids and liquid/cytosol
# Select only cell fractions
plot.data.cell <- plot.data[plot.data$phase %in% c("fract.lip.cell",
                                                   "fract.dis.cell",
                                                   "fract.prot.cell"), ]

# Organize compounds name
plot.data.cell$congener <- factor(plot.data.cell$congener,
                                  levels = c('PCB3', '4-OH-PCB3', '4-OH-PCB3 sulfate',
                                             'PCB11', '4-OH-PCB11', '4-OH-PCB11 sulfate',
                                             'PCB25', '4-OH-PCB25', '4-OH-PCB25 sulfate',
                                             'PCB52', '4-OH-PCB52', '4-OH-PCB52 sulfate'))
# Organize phases
plot.data.cell$phase <- factor(plot.data.cell$phase,
                               levels = c('fract.lip.cell', 'fract.dis.cell',
                                          'fract.prot.cell'))

ggplot(plot.data.cell, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("lipid" , "cytosol", "protein"),
                    values = c("cornsilk3", "brown3", "darkslategray")) +
  theme_classic() +
  theme(aspect.ratio = 10/10,
        text = element_text(size = 14, family = "Helvetica",
                            face = "bold", color = "black")) +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in cell")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
