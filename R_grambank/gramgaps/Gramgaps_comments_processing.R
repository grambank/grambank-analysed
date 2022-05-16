source("requirements.R")

#script written by Alena Witzlack-Makarevich

## Load values.csv
gaps <- read.delim("../grambank/cldf/values.csv", sep = ",")
# length(unique(gaps$Comment)) #42276 unique comments

# Create the copy of the column
gaps <- gaps %>% mutate(Comment_adjusted = Comment)

## basic string cleaning
gaps$Comment_adjusted <- gsub('[[:punct:] ]+',' ',gaps$Comment_adjusted) # remove punctuation

gaps$Comment_adjusted %>% 
  str_to_lower(locale = "en") %>% # convert to low case
  str_trim() %>% # remove whitespace & reduces repeated whitespace
  str_squish() -> gaps$Comment_adjusted # 
# length(unique(gaps$Comment_adjusted)) # 41464

gaps %>% 
  mutate(Comment_class = Comment_adjusted)  -> gaps

# there are comments made up of a number (probably example or table)
# we do not want to remove all numbers in comments, only these ones
gaps %>% 
  mutate(Comment_class = replace(Comment_class, 
                                 which(Comment_class 
                                       %in% c("0", "1", "1 1", "2", "1 3", "10 1 11", "1 16", "1 19", "105", "108 110", "1a e", 
                                              "114", "122 99", "127", "132 135", "145 146", "148", "15", "151 153", 
                                              "158 181", "17", "17 20 21", "18", "18 22", "221", "23", "24", "241 243",
                                              "248 9", "25", "250", "258", "272", "273", "285", "3 1", "3 21", "3 25", "2 447 450",
                                              "212", "259 262", "3 14", "3 4", "3 42", "3 48 3 49", "3 5", "35 39", "3 3 c",
                                              "3 1 b d", "3 25 c", "3 29 a b", "3 3 a", "33 c", "3 49 a", "3 9 b", "31", "338", "35", 
                                              "39", "35", "354 356", "360 368", "38", "39", "39 40", "4", "4 11 6 13", "4 14", "4 5 b c",
                                              "4 14 4 15", "4 33", "4 5", "40", "41", "4 64", "4 63 4 65", 
                                              "412", "44 46",
                                              "440", "45", "46", "47", "47 48", "49", "49 50", "5 10 a", "5 10 b", "52", 
                                              "53 54", "56", "57 58", "6", "6 22 6 23 6 24", "6 59 6 60 6 61", "6 88", "6 89", "7", "7 1", 
                                              "7 14 7 15 7 16", "7 5 2", "7 6", "71 72", "73 74 3", "75", "78", "88", "96", "99 100"
                                       )), 
                                 NA)) -> gaps

# length(unique(gaps$Comment_class)) # 41361

# delete all sorts of frequent irrelevant comments
gaps %>% 
  mutate(Comment_class = replace(Comment_class, 
                                 which(Comment_class 
                                       %in% c("", # there are empty strings
                                              "â", "gr","datapoint inherited from danielsen s work in the database of south american indigenous language structures sails",
                                              "lrc: 206", "false", "lrc 206", 
                                              "initial hg mapping erroneous coding double checked by gb coder",
                                              "probably not", "probably yes", "probably no", 
                                              "revision added", "imported data checked by gb coder manually and value confirmed",
                                              "further confirmation on this datapoint is highly appreciated",
                                              "mm entered source and value", "previously error in import from hg value re evaluated by gb coder")), 
                                 NA)) -> gaps

gaps$Comment_class <- str_replace(gaps$Comment_class, ".*autotranslated.*", replacement = NA_character_)

# length(unique(gaps$Comment_class)) # 40338

# CLASSIFY COMMENTS 
## NOT MENTIONED
gaps$Comment_class <- gaps$Comment_class %>%
  str_replace(".*(not|nor|none|non) (mention|discussed|dicsussed|found|addressed|treated|investigated|covered|described|exhaustive|explicitly|explored|enough).*", "not mentioned") %>%
  str_replace(".*(no|lack of) (description|discussion|example|mention|information|clear|section|list of|evidence).*", "not mentioned") %>% 
  str_replace(".*(insufficient|absent in the grammar|not otherwise mentioned|not overtly discussed|information needed|little information|not well described|short grammar|unable to find|limited data).*", "not mentioned") %>%
  str_replace(".*no .*paradigm.*", "not mentioned") %>%
  str_replace(".*inconclusive.*", "not mentioned") %>%
  str_replace(".*assumed.*absent.*", "not mentioned") %>%
  str_replace(".*only one example.*", "not mentioned") %>%
  str_replace(".*only .* (covered|exemplified).*", "not mentioned") -> gaps$Comment_class

# length(unique(gaps$Comment_class)) #34861

gaps$Comment_class[gaps$Comment_class %in% c("unclear", "not clear", "not quite clear", "not sufficiently discussed",
                                             "not mentionned", #note the typo!
                                             "unmentioned", "non detected", 
                                             "unknown", "not known", 
                                             "no appropriate examples", "not in examples", # TODO is this assignment correct?
                                             "insufficient evidence to say definitively", 
                                             # ranting about the grammars
                                             "the only source is very brief very few example sentences no glossing", 
                                             "restricted material", "not clear and no sentence examples",
                                             "very short unglossed grammar", "very restricted material", 
                                             "one grammar is very short and the other is not glossed with no descriptions unable to determine",
                                             "there is not evidence for or against this in the grammar", 
                                             "pred possession not covered no examples found", "predicate possession not treated", 
                                             "relative clauses are not treated in this source perhaps outside the scope of the descriptive material", 
                                             "no examples found", "probably not no examples found in consulted grammars", 
                                             "needs research", "has to be investigated", 
                                             "data very limited unknown", "unable to tell from the glossing and no overt mention", 
                                             "suppletion not treated in this source perhaps outside the scope of the descriptive material", 
                                             "not evidenced in examples", "not in evidence but difficult to exclude",
                                             "not decribed in the relevant section presumed to be absent", 
                                             "examples do not provide evidence", "only one example found", "not touched upon",
                                             "unable to confirm grammar does not mention interrogation or have relevant examples",
                                             "yes no questions are discussed at the end only the answers to yes no questions are provided no examples given", 
                                             "nominal number marking is not explained", 
                                             "no interrogative information")] <- "not mentioned"

### PASSIM
gaps$Comment_class <- str_replace(gaps$Comment_class, ".*passim.*", "passim")

### NO CATEGORY class 
gaps$Comment_class %>% 
  str_replace(".*(no|any|isn t) (such|morphological case|marking|gender|grammatical gender|noun class|nc gender|verbal person|article|passive|noun class|category|tense|suppletion|comparative construction|grammatical gender|plural marking|class system|dual|inflections|case marking|number marking|person indexing|argument index|index|evidence).*", 
              "no category") %>%
  str_replace(".*not (appear|have).*noun class.*", "no category") %>%
  str_replace(".*no longer functional.*", "no category") %>%
  str_replace(".*not inflected for.*", "no category") %>%
  str_replace(".*there (is|are) no.*", "no category") %>%
  str_replace(".*there does not seem to be.*", "no category") %>%
  str_replace(".*assumed to be absent.*", "no category") -> gaps$Comment_class

#length(unique(gaps$Comment_class))#32735
gaps$Comment_class[gaps$Comment_class %in%
                     c("none found", "none detected", "indication of sex by lexical means is discussed but nothing is said about grammatical gender or noun classes",
                       "nouns do not take any affixes", "no longer functional",
                       "relevant discussion suggests not",
                       "there is no such variation", 
                       "grammatical gender is absent in laboya",
                       "gender is marked on the noun no agreement in the with np dependents",
                       "the only kind of gender that is discussed is lexical gender e g gender of animals as expressed with words like male female",
                       "neither person number tense nor any other grammatical category is coded on the verb these categories may be expressed by independent lexical items if required"
                     )] <- "no category"

#length(unique(gaps$Comment_class))#32727

### SPECIFIC class 
gaps$Comment_class %>% 
  str_replace(".*so called.*", "specific") %>% 
  str_replace("called.*", "specific") %>% 
  str_replace("see e g .*", "specific") %>% 
  str_replace("e g .*", "specific") %>% 
  str_replace(".*unclear (if|whether).*", "specific") %>% 
  str_replace(".*(not|never) (marked|productive).*", "specific") %>% 
  str_replace(".*dedicated number.*", "specific") -> gaps$Comment_class

# length(unique(gaps$Comment_class))#31727

gaps$Comment_class[gaps$Comment_class %in% c("author states there is a possible example of an article but further investigating needs to be done",
                                             "adjectives", "with reflexive", 
                                             "numeration above three is unusual in australian languages capell coate 1984 154 cf also dixon 2004 67",
                                             "zero marking", "zero marked", "zero for inanimates", # check maybe can be grouped under zero
                                             "basic constituent order is sov", "there is only one morphological number marker a freestnading plural particle",
                                             "not productive", "neutral alignment no marking", "core participants are not marked on the verb", 
                                             "there is variation due to topicalisation", 
                                             "only has positive and intensive degree", "subject prefix is the only affix in this language",
                                             "absent from all kimberley languages", "causative",
                                             "not on the verb", "reduplication", "not obligatory", 
                                             "grammar states there is only 9 verbal affixes", 
                                             "note: kim mun share in most part the same grammatical structure and features as biao min [bje] (p338)", 
                                             "this language has extremely free word order",
                                             "no tone", "unclear if productive", 
                                             "tam is indicated through lexical markers", 
                                             "clitics", "clitic", "suffix", "non past", 
                                             "adverbs", "causative prefix", "marked on the verb",
                                             "it doesn t appear so from the listed affixes within the language",
                                             "predicative possession expressed with habeo verb",
                                             "less frequently though not uncommonly words of other parts of speech especially nouns and adverbs but sometimes interjections and particles occur instead of uvs uninflecting verbs in compound verb constructions mcgregor 2004 174",
                                             "most word orders are the same in all bai dialects which is svo note bai has a marked sov word order in a negative or interrogative sentence or when a pronoun is used as one of the two objects",
                                             "the s and a participant is marked by agreement on the predicate in case of a verbal predicate person is expressed by suffixes in the perfect tense and by prefixes and sometimes suffixes in the imperfect tense",
                                             "plural number is only marked in the case markers of proper nouns",
                                             "verbs do not have the category of person and tense", "fingers and toes",
                                             "although the author states that there are three nominal classes body parts kinship terms plus the rest they only differ with respect to their possessive marking this seems to be rather a distinction between alienable and inalienable possession than a distinction between grammatical genders nominal classes in the traditional sense",
                                             "whenever temporal circumstance is expressed this is either done with an elaborate statement of the date or with a phrase like that was the time when"
)] <- "specific"

# length(unique(gaps$Comment_class))#31689

## FROM EXAMPLES ##############
gaps$Comment_class %>%
  str_replace(".*(from|based on|see|various) example.*", "not mentioned") %>%
  str_replace(".*(deduced|inferred) from.*", "not mentioned") %>% 
  str_replace(".*from the example.*", "not mentioned") %>%
  str_replace(".*would have been apparent.*", "not mentioned") -> gaps$Comment_class

#length(unique(gaps$Comment_class))#31345
gaps$Comment_class[gaps$Comment_class %in% c("cf examples", "in examples", "enough examples to exclude", "this is based on one example", 
                                             "would have been apparent in the example sentences", 
                                             "examples support this feature but it is not straightforwardly discussed"
)] <- "not mentioned"

### PAGE MISSING convert to NA
gaps$Comment_class <- str_replace(gaps$Comment_class, ".*page.*missing.*", replacement = NA_character_)
#length(unique(gaps$Comment_class))#31332

### NOTE ON REFERENCES class
gaps$Comment_class %>%
  str_replace(".*also see.*", "note on references or variety") %>%
  str_replace(".*based on.*", "note on references or variety") %>%
  str_replace("see also mc.*", "note on references or variety") %>%
  str_replace(".*see (mc|gb|above|comment|fig).*", "note on references or variety") %>%
  str_replace(".*rumsey.*", "note on references or variety") %>%
  str_replace(".*brownie.*", "note on references or variety") %>%
  str_replace(".*kim mun share.*", "note on references or variety") %>%
  str_replace(".*chapter \\d.*", "note on references or variety") %>%
  str_replace(".*cf gb.*", "note on references or variety") %>%
  str_replace("gr p .*", "note on references or variety") %>%
  str_replace(".*see especially .*", "note on references or variety") %>%
  str_replace(".*table.*", "note on references or variety") %>%
  str_replace("see .* on page.*", "note on references or variety") -> gaps$Comment_class

# length(unique(gaps$Comment_class))#30558

gaps$Comment_class[gaps$Comment_class %in%
                     c("present in all kimberley languages", "lesley stirling pers comm", "examples not glossed"
                     )] <- "note on references or variety"

# length(unique(gaps$Comment_class))#30554

## Classify the remaining long comments as "specific", all the desired classes are shorter than 40 characters
gaps %>% 
  mutate(Comment_length = str_length(Comment_class)) -> gaps

#count the number of characters in the comments and use the information to replace the long ones with "specific"
gaps$Comment_class[gaps$Comment_length > 40] <- "specific"

# length(unique(gaps$Comment_class)) #7862

gaps$Comment_class <- as.factor(gaps$Comment_class)

# convert all infrequent comments to "specific"
gaps %>%
  mutate(Comment_class = fct_lump_min(Comment_class, 11, other_level='specific'))  -> gaps
# length(unique(gaps$Comment_class)) #6

# add a variable where NA is converted to string
gaps %>%
  mutate(Comment_class_NA = str_replace_na(gaps$Comment_class)) -> gaps

gaps$Comment_class_NA <- as.factor(gaps$Comment_class_NA) 

# skip variables we won't use
gaps %>% 
  dplyr::select(ID, Language_ID, Parameter_ID, Value, Source, Comment, Comment_class, Comment_class_NA) -> gaps

OUTPUTDIR <- "output/gramgaps_data/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

# save 
save(gaps, file = paste0(OUTPUTDIR, "/values_tagged_for_comment_class.RData"))
