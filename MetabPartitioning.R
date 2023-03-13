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
logKlip.w <- d$logKlipid.water
logKpro.w <- d$logKprotein.water

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
  V.water.cell <- C.water.cell*Cell/10^6 # L water inside cell
  
  # Temperature correction for Kaw
  R <- 8.3144 # J/(mol K)
  Tst <- 25 # Standard temperature C
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  T2 <- 37 # Incubation temperature C
  T2 <- 273.15 + T2
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/T2-1/Tst.1)) # Ka.w corrected by water and air temps
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculation
  #C.t <- C.dis*(1 + 10^(log.Klip.w)*C.lipid + 10^(log.Kpro.w)*C.prot+
  #                10^(log.Kw.air)*Va/(Vm-V.water.cell)
  den <- 1 + 10^(logKlip.w)*C.lipid + 10^(logKpro.w)*C.prot+
    10^(logKa.w.t)*Va/(Vm-V.water.cell)
  f.dis <- 1/den
  f.dis.m <- f.dis*(Vm-V.water.cell)/Vm
  f.dis.c <- f.dis*V.water.cell/Vm
  f.lip <- 10^(logKlip.w)*C.lipid/den
  f.prot <- 10^(logKpro.w)*C.prot/den
  f.air <- 10^(logKa.w.t)*Va/Vm/den
  f.cell <- f.dis.c + f.lip + f.prot
  f.dis.c.c <- f.dis.c/f.cell
  f.lip.c <- f.lip/f.cell
  f.prot.c <- f.prot/f.cell
  frac <- c(f.dis.m, f.dis.c, f.lip, f.prot, f.air, f.cell, f.dis.c.c,
            f.lip.c, f.prot.c)
}

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener) {
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKlip.w[i], logKpro.w[i]))
}

final.result <- data.frame(congener, result)
names(final.result) <- c("congener", "fract.diss.m", "fract.diss.cell",
                         "fract.lip", "fract.prot", "fract.air",
                         "fract.cell", "fract.diss.cell.cell",
                         "fract.lip.cell", "fract.prot.cell")

# Plots -------------------------------------------------------------------

# (1) Just fractions in cell, medium and air
# create data.frame with needed fractions
p.1 <- final.result[,!names(final.result) %in% c("fract.diss.cell",
                              "fract.lip", "fract.prot",
                              "fract.diss.cell.cell",
                              "fract.lip.cell", "fract.prot.cell")]

# Transform data.frame p.1 to 3 column data.frame
p.1.plot <- melt(p.1, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.1.plot$congener <- factor(p.1.plot$congener,
                           levels = c('PCB3', '4-OH-PCB3', '4-OH-PCB3 sulfate',
                                      'PCB11', '4-OH-PCB11', '4-OH-PCB11 sulfate',
                                      'PCB25', '4-OH-PCB25', '4-OH-PCB25 sulfate',
                                      'PCB52', '4-OH-PCB52', '4-OH-PCB52 sulfate'))

# Organize fraction to be displayed in plot
p.1.plot$phase <- factor(p.1.plot$phase,
                        levels = c('fract.air', 'fract.diss.m',
                                   'fract.cell'))

ggplot(p.1.plot, aes(x = congener, y = fraction, fill = phase)) + 
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
# create data.frame with needed fractions
p.2 <- final.result[,!names(final.result) %in% c("fract.diss.m", "fract.diss.cell",
                                                 "fract.lip", "fract.prot",
                                                 "fract.air", "fract.cell")]
  
# Transform data.frame p.1 to 3 column data.frame
p.2.plot <- melt(p.2, id.var = c("congener"),
                 variable.name = "phase", value.name = "fraction")

# Name the compounds
p.2.plot$congener <- factor(p.2.plot$congener,
                            levels = c('PCB3', '4-OH-PCB3', '4-OH-PCB3 sulfate',
                                       'PCB11', '4-OH-PCB11', '4-OH-PCB11 sulfate',
                                       'PCB25', '4-OH-PCB25', '4-OH-PCB25 sulfate',
                                       'PCB52', '4-OH-PCB52', '4-OH-PCB52 sulfate'))

# Organize fraction to be displayed in plot
p.2.plot$phase <- factor(p.2.plot$phase,
                         levels = c('fract.lip.cell', 'fract.diss.cell.cell',
                                    'fract.prot.cell'))

ggplot(p.2.plot, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("lipid" , "liquid/cytosol", "protein"),
                    values = c("cornsilk3", "brown3", "darkslategray")) +
  theme_classic() +
  theme(aspect.ratio = 10/10,
        text = element_text(size = 14, family = "Helvetica",
                            face = "bold", color = "black")) +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in cell")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
