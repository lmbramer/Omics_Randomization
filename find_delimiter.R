find_delimiter <- function(examp.name){
## this function takes in a character string and returns the most frequent delimeter ##

## store the number of characters in the string ##
num = nchar(examp.name)

## pull each individual character ##
str.parse = substring(examp.name,1:num,1:num)

## identify non-alpha-numeric character locations ##
punct.chars.id = grep("[^[:alnum:]]",str.parse)

## pull non-alpha-numerics ##
punct.chars = str.parse[punct.chars.id]

## count occurences of each character type and return most frequently occuring character ##
seperator = dimnames(table(punct.chars))$punct.chars[which.max(table(punct.chars))]

return(seperator)
}

