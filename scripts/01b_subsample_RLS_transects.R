# remove temperate transects

RLS <- read.csv("data/RLS/RLS_species_NEW_all.csv", sep=";")
RLS<- RLS %>%
  filter(Realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))

# Select 1 year per survey RLS
RLS$Year <- substr(RLS$SurveyDate, 1, 4)

station <- unique(RLS$Station)
RLS_sub <- data.frame()
RLS_st <- vector("list")
for (i in 1:length(station)) {
  RLS_st[[i]] <- RLS %>%
    filter(Station==station[i]) %>%
    filter(Year == max(Year))
  RLS_sub <- rbind(RLS_sub, RLS_st[[i]])
}

RLS_sub <- select(RLS_sub, -Year)  

write.csv(RLS_sub, "data/RLS/RLS_species_NEW.csv")


# filter coral cover
RLS_CC <- read.csv("data/RLS/RLS PQ.csv")

survey <- unique(RLS_CC$SurveyID)
RLS_coral <- data.frame()
Survey_CC <- data.frame(SurveyID=numeric(), coral_cover=numeric())
for (i in 1:length(survey)) {
  RLS_coral <- RLS_CC %>%
    filter(SurveyID==survey[i])%>%
    filter(MajorCategory=="Coral")
  Survey_CC[i,1] <- survey[i]
  Survey_CC[i,2] <- sum(RLS_coral$Percent.Coverage)
}

RLS_sub_CC <- left_join(RLS_sub, Survey_CC, by="SurveyID")

write.csv(RLS_sub_CC, "data/RLS/RLS_species_NEW.csv")
