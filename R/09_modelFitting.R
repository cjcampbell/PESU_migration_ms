
# Setup -------------------------------------------------------------------

library(MuMIn)
library(car)
options(na.action = "na.fail")

myResults <- read.csv( file.path(wd$bin, "myResults.csv") )

df <- myResults %>% 
  filter(Season == "Winter") %>% 
  filter(!is.na(Sex), !is.na(Size), !is.na(Region_long)) %>% 
  dplyr::mutate(Region_long = factor(Region_long, levels = rev(c("North-central", "Northwest"))))


# Fit global model --------------------------------------------------------

# Check out predictor formats and proportions.
df$Sex  %>% as.character %>% as.factor %>% as.numeric %>% hist
df$Size %>% as.character %>% as.factor %>% as.numeric %>% hist
df$Region_long %>% as.character %>% as.factor %>% as.numeric %>% hist

# Check out response dist.
df$probDistTraveledAtThreshold_0.25 %>% hist

# Fit model.
m1 <- glm(probDistTraveledAtThreshold_0.25~Sex+Size+Region_long, data = df)
d1 <- dredge(m1)
m2 <- glm(probDistTraveledAtThreshold_0.25~Size+Region_long, data = df)
d2 <- dredge(m2)
car::vif(m2)


# Put model estimates into temporary data.frames:
model1Frame <- data.frame(Variable = rownames(summary(m2)$coef),
                          Coefficient = summary(m2)$coef[, 1],
                          SE = summary(m2)$coef[, 2],
                          modelName = "Top Model", row.names = NULL)
model2Frame <- data.frame(Variable = rownames(summary(m1)$coef),
                          Coefficient = summary(m1)$coef[, 1],
                          SE = summary(m1)$coef[, 2],
                          modelName = "Global Model",  row.names = NULL)

# Combine these data.frames
allModeldf <- data.frame(rbind(model1Frame, model2Frame)) %>% 
  dplyr::mutate(
    modelName = factor(modelName, levels = c("Global Model","Top Model")),
    Variable = case_when(
      Variable == "Region_longNorth-central" ~ "Karst Region (N-C)",
      Variable == "SizeSmall" ~ "Colony Size (Small)",
      Variable == "SexMale" ~ "Sex (Male)",
      TRUE ~ Variable
      ),
      Variable = factor(Variable, levels = rev(c("Karst Region (N-C)", "Colony Size (Small)", "Sex (Male)", "(Intercept)")))
    )

save(m2, m1, allModeldf,  file = file.path(wd$bin, "modelResults.Rdata"))


# check model results  -----------------------------------------------------------
library(gtsummary)
library(performance)

# How much variance is explained?
performance::r2(m2) 

# Table
gtsummary::tbl_regression(m2)

# Plot resids.
plot(m2)
sjPlot::plot_model(m2)
