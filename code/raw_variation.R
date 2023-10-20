################################################################################
# Variation in blood feeding behavior
################################################################################


# Load packages and functions ---------------------------------------------
library(tidyverse)


# Read in raw data --------------------------------------------------------
filenames <- list.files("./data/raw/", pattern = "*.csv", 
                        full.names = TRUE, recursive = TRUE)


ldf <- lapply(filenames, read.csv)

data <- tibble(T = double(), biting_rate = double(), ref = character(), 
               mosquito_species = character(), parasite = character(),
               dataset_ID = integer())

for (i in 1:length(ldf)) {
  
  cols_to_rename <- c(trait.name = "Trait.Name", 
                      mosquito_species = "mosquito.species",
                      mosquito_species = "host.code",
                      mosquito_species = "specie",
                      parasite = "pathogen", 
                      parasite = "paras.code",
                      ref = "Ref"
                      )
  
  temp_frame <- ldf[[i]] %>% 
    rename(any_of(cols_to_rename)) %>% 
    select(any_of(c("trait.name", "T", "trait", "ref", "mosquito_species", "parasite"))) %>% 
    # assign parasite and mosquito species for data sets where they aren't included in the actual file
    # needs to be done for datasets: 3, 4, 
    mutate(dataset_ID = i,
      parasite = ifelse("parasite" %in% names(.), parasite, NA),
      mosquito_species = ifelse("mosquito_species" %in% names(.), mosquito_species, NA),) %>% 
    filter(trait.name %in% c("a", "GCR", "GCD", "gcr", "gcd")) %>%
    mutate(biting_rate = case_match(trait.name,
                                    "a" ~ trait,
                                    "GCR" ~ trait,
                                    "gcr" ~ trait,
                                    "GCD" ~ 1/trait,
                                    "gcd" ~ 1/trait
    )
    ) %>% 
    select(-c(trait.name, trait))
  
  data <- rbind(data, temp_frame)
}

res <- lapply(ldf, summary)
names(res) <- substr(filenames, 6, 30)


# Visualize variation -----------------------------------------------------


