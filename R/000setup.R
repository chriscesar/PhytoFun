# 000setup.R ####
# Set up project themes/metadata

# Colours ####

cbPalette <- c( ### colourblind-friendly chart colour palette
  # comments reflect Red/Green/Blue equivalents
  "#0072B2", #000, 114, 178
  "#e79f00", #231, 159, 0
  "#009E73", #000, 158, 115
  "#9ad0f3", #154, 208, 243
  "#000000", #0, 0, 0
  "#D55E00", #213, 94, 0
  "#CC79A7", #204, 121, 167
  "#DEAAC6", #222, 170, 198 - The Wash
  "#F0E442"  #240, 228, 66
)
cbPaletteTxt <- c(
  "#0072B2", #000, 114, 178
  "#e79f00", #231, 159, 0
  "#009E73", #000, 158, 115
  "#78b0d1", #120, 176, 209
  "#000000", #0, 0, 0
  "#D55E00", #213, 94, 0
  "#CC79A7", #204, 121, 167
  "#F0E442") #240, 228, 66

cbPaletteshr <- c("#F15854","#B276B2","#B2912F","#4D4D4D")

cbPaletteFill <- c( ### colourblind-friendly chart colour palette - lighter cols for fills
  # comments reflect Red/Green/Blue equivalents
  "#238ACD", #000, 114, 178
  "#FFB935", #231, 159, 0
  "#23B888", #000, 158, 115
  "#B5EAFE", #154, 208, 243
  "#232323", #0, 0, 0
  "#F07823", #213, 94, 0
  "#E894BC", #204, 121, 167
  "#F8C4E1", #222, 170, 198 - The Wash
  "#FFFE77"  #240, 228, 66
)

cbPalette2 <- c("#999999",
                "#E69F00",
                "#56B4E9",
                "#CC79A7",
                "#009E73",
                "#F0E442",
                "#0072B2",
                "#D55E00",
                "#444444",
                "#C34D55",
                "#33A2C4",
                "#554C31",
                "#C5C221",
                "#5531A1",
                "#B32C55",
                "#BB3593")

ppi <- 300 #figure resolution

# Plot themes ####
ggplot2::theme_set(ggthemes::theme_few())
theme_use <- theme(
  strip.text = element_text(face = 2),
  axis.title.x = element_blank(),
  axis.title.y = element_text(face=2),
  axis.text = element_text(face=2),
  axis.text.x = element_text(size = 12),
  plot.caption = ggtext::element_markdown(face=2,size=12),
  plot.title = element_text(face=2,size = 14),
  
)

# theme( # use theme_get() to see available options
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_blank(),
#     legend.key = element_blank(), # switch off the rectangle around symbols in the legend
#     # legend.position = ("bottom"),
#     legend.justification = (NULL),
#     # legend.direction=("horizontal"),
#     legend.position="topright",
#     legend.text = element_text(size=16,family="TN"),
#     legend.title = element_text(size=16,family="TN"),
#     text=element_text(size=16,family="TN") ,
#     strip.text.y = element_text(size = 13,angle = 0,family="TN"),
#     axis.title.x=element_text(vjust=0.3, size=16,family="TN"),
#     axis.title.y=element_text(vjust=0.6, angle=90, size=16,family="TN"),
#     axis.text.x=element_text(size=16,family="TN"),
#     axis.text.y=element_text(size=16,family="TN"),
#     axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
#     axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
#     strip.background = element_blank())
# 
#           

# Folders ####
source("R/dataFolders.R")