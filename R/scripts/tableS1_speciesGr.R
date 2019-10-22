#### TABLE SPECIES GROUP #####

### PACKAGES ####

library(dplyr)
library(knitr)

### DATA ####

tree_code <- read.csv2("data/ref_spCode.csv")

### GROUP ####

pioneer <- c("BETPAP", "BETPOP", 
             "POPBAL", "POPDEL", "POPGRA", "POPTRE", 
             "PRUPEN", "SALSP", "SORSP")

temperate <- c("ACENEG", "ACENIG", "ACEPEN", "ACERIN", "ACERUB", "ACESAC", "ACESPI", 
               "AMESP", "BETALL",
               "CARCAR", "CARCOR", "CAROVA", "FAGGRA", 
               "FRAAME", "FRANIG", "FRAPEN", 
               "JUGCIN", "OSTVIR",
               "PICRUB", "PINRES", "PINRIG", "PINSTR",
               "PRUSER", 
               "QUEALB", "QUEBIC", "QUEMAC", "QUERUB", 
               "THUOCC", "TILAME", "TSUCAN", 
               "ULMAME", "ULMRUB", "ULMTHO")

boreal <- c("ABIBAL", "LARLAR", "PICGLA", "PICMAR", "PINBAN")


tree_code <- tree_code %>% filter(spCode %in% c(temperate,boreal,pioneer)) %>%
  mutate(Group = case_when(spCode %in% pioneer ~ "Pioneer",
                           spCode %in% boreal ~ "Boreal",
                           spCode %in% temperate ~ "Temperate")) %>% 
  arrange(Group, complete.name)

tree_code <- tree_code %>%
  arrange(Group) %>%
  dplyr::select(complete.name, VERNACULAR) %>%
  rename("Species name" = complete.name, "Vernacular name" = VERNACULAR)


kable(tree_code, format = 'latex', longtable=T, booktabs = T, linesep = "", escape = FALSE) %>%
  column_spec(1, italic = T) %>%
  row_spec(0, bold = TRUE) %>%
  group_rows("Temperate", 15, 47) %>% 
  group_rows("Boreal", 1, 5) %>%
  group_rows("Pioneer", 6,14) %>%
  kable_styling(latex_options =c("repeat_header"))

