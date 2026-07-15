# utility functions for cleaning VCCC data

# classify ART combination regimen into NRTI, PI, NNRTI, INSTI, Other
classify_art_regimen <- function(regimen) {
  case_when(
    grepl("lamivudine", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("zidovudine", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("nelfinavir", regimen, ignore.case = TRUE) ~ "PI",
    grepl("DIDANOSINE", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("EFAVIRENZ", regimen, ignore.case = TRUE) ~ "NNRTI",
    grepl("ABACAVIR", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("INDINAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("ATAZANAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("RITONAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("TENOFOVIR", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("EMTRICITABINE", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("STAVUDINE", regimen, ignore.case = TRUE) ~ "NRTI",
    grepl("NEVIRAPINE", regimen, ignore.case = TRUE) ~ "NNRTI",
    grepl("RALTEGRAVIR", regimen, ignore.case = TRUE) ~ "INSTI",
    grepl("BLINDED ART", regimen, ignore.case = TRUE) ~ "OTHER",
    grepl("DARUNAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("LOPINAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("RILPIVIRINE", regimen, ignore.case = TRUE) ~ "NNRTI",
    grepl("FOSAMPRENAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("ENFUVIRTIDE", regimen, ignore.case = TRUE) ~ "OTHER",
    grepl("SAQUINAVIR", regimen, ignore.case = TRUE) ~ "PI",
    grepl("ETRAVIRINE", regimen, ignore.case = TRUE) ~ "NNRTI",
    grepl("Hx ARV Treatment, Med unknown", regimen, ignore.case = TRUE) ~ "OTHER",
    grepl("MARAVIROC", regimen, ignore.case = TRUE) ~ "OTHER",
    grepl("DOLUTEGRAVIR", regimen, ignore.case = TRUE) ~ "INSTI",
    .default = NA
  )
}

# classify ART regimen into NNRTI-containing vs not
regimen_contains_nnrti <- function(regimen) {
  as.numeric(
    grepl("EFAVIRENZ", regimen, ignore.case = TRUE) | grepl("NEVIRAPINE", regimen, ignore.case = TRUE) |
      grepl("RILPIVIRINE", regimen, ignore.case = TRUE) | grepl("ETRAVIRINE", regimen, ignore.case = TRUE)
  )
}
regimen_contains_nrti <- function(regimen) {
  as.numeric(
    grepl("lamivudine", regimen, ignore.case = TRUE) | grepl("zidovudine", regimen, ignore.case = TRUE) |
    grepl("DIDANOSINE", regimen, ignore.case = TRUE) | grepl("ABACAVIR", regimen, ignore.case = TRUE) |
    grepl("TENOFOVIR", regimen, ignore.case = TRUE) | grepl("EMTRICITABINE", regimen, ignore.case = TRUE) |
    grepl("STAVUDINE", regimen, ignore.case = TRUE)
  )
}
regimen_contains_pi <- function(regimen) {
  as.numeric(
    grepl("nelfinavir", regimen, ignore.case = TRUE) | grepl("INDINAVIR", regimen, ignore.case = TRUE) |
      grepl("ATAZANAVIR", regimen, ignore.case = TRUE) | grepl("RITONAVIR", regimen, ignore.case = TRUE) |
      grepl("DARUNAVIR", regimen, ignore.case = TRUE) | grepl("LOPINAVIR", regimen, ignore.case = TRUE) |
      grepl("FOSAMPRENAVIR", regimen, ignore.case = TRUE) | grepl("SAQUINAVIR", regimen, ignore.case = TRUE)
  )
}
regimen_contains_insti <- function(regimen) {
  as.numeric(
    grepl("RALTEGRAVIR", regimen, ignore.case = TRUE) | grepl("DOLUTEGRAVIR", regimen, ignore.case = TRUE)
  )
}


# get first nonvalid and valid ART start date and CD4 count after enrollment
getFirstARTPostEnrollment <- function (sub_id, art_df, demo_df) {

  art_df <- art_df %>% dplyr::filter(CFAR_PID == sub_id) %>%
    dplyr::arrange(AGE_AT_MED_START)
  demo_df <- demo_df %>% dplyr::filter(CFAR_PID == sub_id)

  firstOccur <- art_df$AGE_AT_MED_START >= demo_df$AGE_AT_FIRST_VISIT

  if (any(firstOccur == FALSE) == TRUE) {
    # if any are FALSE, ie. ART before enrollment, exclude subject
    return(NA)
  } else {
    # find position of FIRST TRUE such that firstOccur is TRUE
    return(art_df[Position(isTRUE, firstOccur),])
  }

}

# return the closest CD4 count and viral load
getClosestCD4toART <- function(sub_id, lab_df, art_df, value = "CD4 COUNT") {

  lab_df <- lab_df %>% dplyr::filter(CFAR_PID == sub_id)
  art_df <- art_df %>% dplyr::filter(CFAR_PID == sub_id)

  lastLab <- (art_df$AGE_AT_MED_START - 180/365 < lab_df$AGE_AT_RESULT_DATE) & (lab_df$AGE_AT_RESULT_DATE <  art_df$AGE_AT_MED_START + 30/365)

  if (sum(lastLab) == 0) {
    return(tibble(AGE_AT_RESULT_DATE = NA, CFAR_PID = sub_id, YEAR_OF_LAB_RESULT = NA, cd4 = NA, cd4_percent = NA, vl = NA))
  } else {
    lab_closeCD4 <- lab_df[lastLab,]
    lab_closeCD4$dist <- lab_closeCD4$AGE_AT_RESULT_DATE - art_df$AGE_AT_MED_START
    lab_closeCD4_wide <- lab_closeCD4 %>%
      mutate(testName = case_when(
        testName == "CD4 COUNT" ~ "cd4",
        testName == "CD4 PERCENT" ~ "cd4_percent",
        .default = "vl"
      )) %>%
      select(-units) %>%
      pivot_wider(names_from = testName, values_from = RESULT_NUMERIC)
    return(lab_closeCD4_wide[which.min(abs(lab_closeCD4_wide$dist)),] %>% select(-dist))
  }

}


getFirstADE <- function(sub_id, demo_df, art_df, lab_df, ade_df) {

  demo_sub <- demo_df %>% dplyr::filter(CFAR_PID == sub_id)
  art_sub <- art_df %>% dplyr::filter(CFAR_PID == sub_id)
  lab_sub <- lab_df %>% dplyr::filter(CFAR_PID == sub_id)

  demo_sub$delta <- NA
  lab_sub$delta <- NA

  if (sub_id %in% unique(ade_df$CFAR_PID) == FALSE) {
    # add last known date of contact if subject does not have opportunistic infection and has follow-up after ART start date
    if (is.na(demo_sub$AGE_AT_DEATH) == TRUE) {
      if (lab_sub[nrow(lab_sub),]$AGE_AT_RESULT_DATE > art_sub$AGE_AT_MED_START) {
        lab_sub[nrow(lab_sub),]$delta <- 0
        lab_sub <- lab_sub %>% dplyr::select(AGE = AGE_AT_RESULT_DATE, CFAR_PID, delta)
        return(lab_sub[nrow(lab_sub),])
      } else {
        return(NA)
      }
    } else {
      demo_sub$delta <- 0
      demo_sub <- demo_sub %>% dplyr::select(AGE = AGE_AT_DEATH, CFAR_PID, delta)
      return(demo_sub)
    }
  } else {

    ade_sub <- ade_df %>% dplyr::filter(CFAR_PID == sub_id)

    ade_sub$delta <- NA
    ade_sub <- ade_sub %>% dplyr::select(AGE = AGE_AT_DX_ONSET, CFAR_PID, delta)

    for(l in 1:nrow(ade_sub)){
      #ignore if ADE was before or at enrollment
      if (ade_sub[l,]$AGE <= demo_sub$AGE_AT_FIRST_VISIT) {
        return(NA)
      }
      # ignore if ADE was before or at ART start
      if (ade_sub[l,]$AGE <= art_sub$AGE_AT_MED_START) {
        return(NA)
      }
      # desired event of OI at least 6 months after ART start
      if (ade_sub[l,]$AGE > 0.5 + art_sub$AGE_AT_MED_START) {
        ade_sub[l,]$delta <- 1
        return(ade_sub[l,])
      }
      # desired event of OI within 6 months of ART start
      if (ade_sub[l,]$AGE <= 0.5 + art_sub$AGE_AT_MED_START) {
        # see if there is another ADE after 6 months of ART start
        for (c in l:nrow(ade_sub)) {
          if (ade_sub[c,]$AGE > 0.5 + art_sub$AGE_AT_MED_START) {
            ade_sub[c,]$delta <- 1
            return(ade_sub[c,])
          }
        }
        # if no other ADE, censor subject
        if (is.na(demo_sub$AGE_AT_DEATH) == TRUE) {
          # if age at death is missing, find if last test date is after ART start and add it as censoring date
          if(ade_sub[nrow(ade_sub),]$AGE > art_sub$AGE_AT_MED_START){
            ade_sub[nrow(ade_sub),]$delta <- 0
            return(ade_sub[nrow(ade_sub),])
          } else {
            return(NA)
          }
        } else {
          # if age at death available, add that as censoring date
          demo_sub$delta <- 0
          demo_sub <- demo_sub %>% dplyr::select(AGE = AGE_AT_DEATH, CFAR_PID, delta)
          return(demo_sub)
        }

      }
    }
  }
}
