################################################################################
# Variation in blood feeding behavior
################################################################################


# Load packages and functions ---------------------------------------------
library(tidyverse)
library(MASS)

# Read in raw data --------------------------------------------------------
filenames <- list.files("./data/raw/biting_rates/", pattern = "*.csv", 
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
    dplyr::select(any_of(c("trait.name", "T", "trait", "ref", "mosquito_species", "parasite"))) %>%
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
    dplyr::select(-c(trait.name, trait))
  
  data <- rbind(data, temp_frame)
}

# Visualize variation -----------------------------------------------------

bad_plot <- data %>% 
  arrange(T) %>% 
  # filter(biting_rate > 1/1000) %>% 
  ggplot(aes(x = T, y = biting_rate, color = mosquito_species)) +
  geom_point() +
  geom_line()



# Scott et al 2000 data ---------------------------------------------------
mbm_freqs <- read.csv("./data/raw/Scott_2000.csv") %>% 
  group_by(Site) %>% 
  mutate(props = Count / sum(Count)) %>% 
  ungroup()

fit_mbm_freqs <- mbm_freqs %>% 
  filter(Site == "Thailand") %>% 
  rename(x = "N_bloodmeals", y = "Count") %>% 
  dplyr::select(-Site)
  
x = c(rep(1,fit_mbm_freqs$y[1]),
      rep(2,fit_mbm_freqs$y[2]),
      rep(3,fit_mbm_freqs$y[3])
      )
exp_fit_Thailand <- fitdistr(x, densfun = "exponential")

exp_curve <- tibble(x = seq(0.0, 4.0, length.out = 100)) %>% 
  mutate(y = exp_fit_Thailand$estimate*exp(-exp_fit_Thailand$estimate * x),
         Site = "Thailand")

mbm_plot <- mbm_freqs %>% 
  ggplot(aes(x = N_bloodmeals, y = Count, fill = Site)) +
  geom_col(position = "dodge") +
  geom_line(data = exp_curve, aes(x = x, y = y*1300, lty = Site))
mbm_plot

