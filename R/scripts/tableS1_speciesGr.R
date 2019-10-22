#### TABLE SPECIES GROUP #####

### PACKAGES ####

library(dplyr)
library(knitr)

### DATA ####

tree_code <- read.csv2("data/ref_spCode.csv")

### GROUP ####

pioneer <- c("BETPAP", "BETPOP", "CRASP",
             "POPBAL", "POPDEL", "POPGRA", "POPTRE", "PRUPEN", "SALSP", "SORSP")

temperate <- c("ACENEG", "ACENIG", "ACEPEN", "ACERIN", "ACERUB", "ACESAC", "ACESPI", 
               "AMESP", "BETALL",
               "CARCAR", "CARCOR", "FAGGRA", "FRAAME", "FRANIG", "FRAPEN", "JUGCIN",
               "OSTVIR",
               "PICRUB", "PINRES", "PINRIG", "PINSTR",
               "PRUSER", "PRUVIR",
               "QUEALB", "QUEBIC", "QUEMAC", "QUERUB", "THUOCC",
               "TILAME", "TSUCAN", "ULMAME", "ULMRUB", "ULMTHO")

boreal <- c("ABIBAL","ALNRUG", "LARLAR", "PICGLA", "PICMAR", "PINBAN")

tree_code <- tree_code %>% filter(spCode %in% c(temperate,boreal,pioneer)) %>%
  mutate(Group = case_when(spCode %in% pioneer ~ "Pioneer",
                           spCode %in% boreal ~ "Boreal",
                           spCode %in% temperate ~ "Temperate"))

tree_code <- tree_code %>%
  arrange(Group) %>%
  dplyr::select(complete.name, VERNACULAR) %>%
  rename("Species name" = complete.name, "Vernacular name" = VERNACULAR)


kable(tree_code, format = 'latex', longtable=T, booktabs = T, linesep = "", escape = FALSE) %>%
  column_spec(1, italic = T) %>%
  row_spec(0, bold = TRUE) %>%
  group_rows("Temperate", 17, 49) %>% 
  group_rows("Boreal", 1, 6) %>%
  group_rows("Pioneer", 7,16) %>%
  kable_styling(latex_options =c("repeat_header"))
  
#%>%
  
  cat(file = "res/TableS1.md", sep = "/n")
