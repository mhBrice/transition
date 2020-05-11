#### TABLE S1. SPECIES GROUP #####

### PACKAGES ####

source("R/functions/packages.R")

### DATA ####

source("R/functions/prep_data.R")

tree_code <- read.csv2("data/ref_spCode.csv")

sp_ba <- readRDS("data/sp_mat_ba_nov2019.RDS") %>% filter(ID_PE %in% states_ba$ID_PE)

### GROUP ####

pioneer <- c("BETPAP", "BETPOP",
             "POPBAL", "POPDEL", "POPGRA", "POPTRE",
             "PRUPEN", "SALSP", "SORSP")

temperate <- c("ACENEG", "ACENIG", "ACEPEN", "ACERIN", "ACERUB", "ACESAC", "ACESPI",
               "AMESP", "BETALL",
               "CARCAR", "CARCOR", "FAGGRA",
               "FRAAME", "FRANIG", "FRAPEN",
               "JUGCIN", "OSTVIR",
               "PICRUB", "PINRES", "PINRIG", "PINSTR",
               "PRUSER",
               "QUEALB", "QUEBIC", "QUEMAC", "QUERUB",
               "THUOCC", "TILAME", "TSUCAN",
               "ULMAME", "ULMRUB", "ULMTHO")

boreal <- c("ABIBAL", "LARLAR", "PICGLA", "PICMAR", "PINBAN")

sp_ba <- sp_ba %>%
  mutate(SORSP = SORAME + SORDEC) %>%
  select(-MALSP, -SORAME, -SORDEC, -ALNCRI, -CRASP, -ALNRUG, -PRUVIR)

Frequency <- sp_ba %>% group_by(ID_PE) %>%
  summarise_at(vars(boreal, pioneer, temperate), function(x) any(x>0))
Frequency <- colSums(Frequency[,-1])
Frequency <- as.data.frame(Frequency)
Frequency$spCode <- row.names(Frequency)

tree_code <- tree_code %>% filter(spCode %in% c(temperate,boreal,pioneer)) %>%
  mutate(Group = case_when(spCode %in% pioneer ~ "Pioneer",
                           spCode %in% boreal ~ "Boreal",
                           spCode %in% temperate ~ "Temperate")) %>%
  arrange(Group, complete.name) %>%
  left_join(Frequency, by = "spCode")

tree_code <- tree_code %>%
  arrange(Group) %>%
  dplyr::select(complete.name, VERNACULAR, Frequency) %>%
  rename("Species name" = complete.name, "Vernacular name" = VERNACULAR)


kable(tree_code, format = 'latex', longtable = TRUE, booktabs = TRUE,
  linesep = "", escape = FALSE) %>%
  column_spec(1, italic = TRUE) %>%
  row_spec(0, bold = TRUE) %>%
  group_rows("Temperate", 15, 46) %>%
  group_rows("Boreal", 1, 5) %>%
  group_rows("Pioneer", 6,14) %>%
  kable_styling(latex_options =c("repeat_header"))

