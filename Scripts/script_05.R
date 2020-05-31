`+`(1, 2, 3)
1 + 2 + 3 
`+`(1, `+`(2, 3))
`%cup%` <- function(x, y) union(x, y)
1:10 %cup% 2:15
1:10 %cup% 2:15 %cup% 16:18
union(1:10, union(2:15, 16:18))
max(1:10)
sprintf("%d", max(1:10))
1:10 %>% { sprintf("%d", .) }
1:10 %>% max %>% { sprintf("%d", .) }
. <- 10
.
head(df(1:1000))
head(data.frame(1:1000))
data.frame(1:1000) %>% head
data.frame(1:1000) 
data.frame(1:1000) %>% head
head(data.frame(1:1000) )
1:1000 %>% data.frame %>% head
df
df[,1]
df %>% `[`
library(tidyverse)
df %>% `[`
df %>% extract(1)
?extract
library(magrittr)
df %>% extract(1)
list(a=1, b=2) %$% print(a)
a
data.frame(id=1:10, cont=LETTERS[1:10])
data.frame(id=1:10, cont=LETTERS[1:10])[,1:2]
data.frame(id=1:10, cont=LETTERS[1:10])[,1]
df <- data.frame(id=1:10, cont=LETTERS[1:10])
sel <- columns(df) %in% c("id", "cont", "lfc")
sel <- colnames(df) %in% c("id", "cont", "lfc")
sel
df[ , sel ]
df <- data.frame(id=1:10, x=LETTERS[1:10])
sel <- colnames(df) %in% c("id", "cont", "lfc")
sel
df[ , sel ]
data.frame(ID=c(1:500, "lost"), stringsAsFactors=T)
data.frame(ID=c(1:500, "lost"), stringsAsFactors=T) %>% head
class(data.frame(ID=c(1:500, "lost"), stringsAsFactors=T)[,1])
data.frame(ID=c(sample(1:500), "lost"), stringsAsFactors=T) %>% head
data.frame(ID=c(sample(1:500), "lost"), stringsAsFactors=T)[,1] * 2
table1
as_tibble(data.frame(1:1000))
df <- data.frame(1:10)
rownames(df) <- LETTERS[1:10]
df
as_tibble(df)
df %>% rownames_to_column(ID) %>% as_tibble
df %>% rownames_to_column("ID") %>% as_tibble
df[ "A", ]
df[ "B", ]
df %>% rownames_to_column("ID") %>% as_tibble %>% filter(ID == "A")
df %>% rownames_to_column("ID") %>% as_tibble %>% filter(ID == "A") %>% pull
df %>% rownames_to_column("ID") %>% as_tibble %>% filter(ID == "A") %>% select("X1.10")
tf <- df %>% rownames_to_column("ID") %>% as_tibble 
tf
fg[,1]
tf[,1]
tf[,1] %>% pull(1)
tf[,1] %>% pull
starwars
print(starwars, n=Inf, width=Inf)
starwars %>% data.frame
starwars[ starwars$films == "Attack of the Clones", ]
starwars %>% filter(films == "Atack of the Clones")
starwars %>% filter(films == "Atack of the Clones") %>% pull(name)
starwars %>% filter(films == "Atack of the Clones") %>% pull("name")
starwars %>% filter(films == "Attack of the Clones") %>% pull(name)
starwars[ starwars$films == "Attack of the Clones", "name" ]

starwars %>% filter(films == "Atack of the Clones") %>% select(name)
starwars %>% filter(films == "Attack of the Clones") %>% select(name)
s.df <- as.data.frame(starwars)
s.df[ , "name", drop=FALSE]
starwars 
starwars %>% filter(films != "Attack of the Clones" & homeworld == "Tatooine") %>% select(name)
starwars %>% filter(films != "Attack of the Clones" & eye_color == "red") %>% select(name)
starwars %>% filter(films != "Attack of the Clones" & eye_color == "red") 
starwars %>% filter(films != "Attack of the Clones" & eye_color == "red") %>% arrange(mass)
starwars %>% filter(films != "Attack of the Clones" & eye_color == "red") %>% arrange(-mass)
starwars %>% filter(films != "Attack of the Clones" & eye_color == "red") %>% arrange(desc(mass))
starwars %>% filter(films == "Attack of the Clones" & eye_color == "red") %>% arrange(desc(mass) )
starwars %>% filter(films == "Attack of the Clones") %>% arrange(desc(mass) )
starwars %>% filter(films == "Attack of the Clones") %>% arrange(homeworld, desc(mass))
starwars %>% filter(films == "Attack of the Clones") %>% arrange(homeworld, desc(height))
starwars %>% filter(films == "Attack of the Clones") %>% arrange(homeworld, gender, eye_color)
starwars %>% filter(films == "Attack of the Clones") %>% arrange(homeworld, gender, eye_color)
sw1 <- starwars[ order(starwars$eye_color), ]
sw2 <- sw2[ order(sw1$gender), ]
sw2 <- sw1[ order(sw1$gender), ]
sw2[ order(sw2$homeworld), ]
sw2 <- sw2[ sw2$films == "Attack of the Clones", ]
sw2[ order(sw2$homeworld), ]
starwars %>% filter(films == "Attack of the Clones") %>% arrange(homeworld, gender, eye_color)
starwars %>% summarise(hair_color)
starwars %>% group_by(hair_color) %>% summarise
starwars %>% group_by(hair_color) 
?summarise
starwars %>% summarise(m=mean(mass, na.rm=TRUE))
starwars %>% summarise(m=mean(mass, na.rm=TRUE), m.var=var(mass, na.rm=TRUE))
starwars %>% arrange(-mass)
starwars %>% summarise(m=mean(mass), m.var=var(mass))
mean(c(1, 2, 3, NA))
mean(c(1, 2, 3, NA), na.rm=TRUE)
starwars %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE))
starwars %>% group_by(species) 
identifal(as.data.frame(starwars), as.data.frame(starwars %>% group_by(species)))
identical(as.data.frame(starwars), as.data.frame(starwars %>% group_by(species)))
starwars %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE))
starwars %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n())
n() 
starwars %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n()) %>% filter(n > 1)
starwars %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n()) %>% filter(n > 1) %>% arrange(-desc(n))

starwars %>% group_by(species) %>% add_tally()
starwars %>% group_by(species) %>% add_count()
?add_tally
starwars %>% add_tally(species)
starwars %>% add_tally("species")
starwars %>% add_tally(~ species)
?add_tally
starwars %>% add_count(species)
mtcars
mtcars %>% select(mpg)
mtcars %>% select(MilesPerGallon=mpg)
mtcars %>% head
mtcars %>% select(mpg:disp)
mtcars %>% select(mpg:hp)
mtcars %>% select(mpg:hp) %>% group_by(cyl) %>% summarise(av_mgp=mean(mpg))

mtcars %>% select(starts_with("c"))
mtcars %>% select(mpg)
mpg <- "cyl"
mtcars %>% select(mpg)
mtcars %>% select(!!mpg)
mtcars %>% mutate(lp100km=282.5/mpg)
?starts_with
?matches
?matches
starwars %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n()) %>% filter(n > 1) %>% arrange(-desc(n))
df <- starwars
starwars %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n()) %>% filter(n > 1) %>% arrange(-desc(n))
df
df %>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n()) %>% filter(n > 1) %>% arrange(-desc(n))
df %<>% group_by(species) %>% summarise(m=mean(mass, na.rm=TRUE), h=mean(height, na.rm=TRUE), n=n()) %>% filter(n > 1) %>% arrange(-desc(n))
df
lapply(tmp, function(x) paste(x, collapse=", "))
lapply(tmp, function(.) paste(., collapse=", "))
map(tmp, ~ paste(., collapse=", "))
tmp
lapply(names(tmp), function(n) paste("The element", n, "has", length(tmp[[n]]), "elements")
lapply(names(tmp), function(n) paste("The element", n, "has", length(tmp[[n]]), "elements") })
lapply(names(tmp), function(n) paste("The element", n, "has", length(tmp[[n]]), "elements"))
lapply(names(tmp), function(n) paste("The element", n, "has", length(tmp[[n]]), "element(s)"))
lapply(tmp, function(n) paste("The element has", length(n), "element(s), but I have no idea what its name is"))
imap(tmp, ~ { paste("The element", .y, "has", length(.x), "element(s)" })
imap(tmp, ~ { paste("The element", .y, "has", length(.x), "element(s)") })
starwars
starwars$vehicles[1]
starwars %>% mutate(n_vehicle=map(vehicles, length))
starwars %>% mutate(n_vehicle=map_int(vehicles, length))
map2_chr(starwars$name, starwars$vehicles, ~ { paste(.x, length(.y)) })
starwars %>% select(name, hair_color, eye_color)
starwars %>% select(name, hair_color, eye_color) %>% mutate_at(~ gsub("none", NA, .))
starwars %>% select(name, hair_color, eye_color) %>% mutate_at(vars(~ gsub("none", NA, .)))
starwars %>% select(name, hair_color, eye_color) %>% mutate_all(~ gsub("none", NA, .))
starwars %>% select(name, hair_color, eye_color) %>% mutate_at(hair_color:eye_color, ~ gsub("none", NA, .))
starwars %>% select(name, hair_color, eye_color) %>% mutate_at(vars(hair_color:eye_color), ~ gsub("none", NA, .))
starwars %>% select(name, hair_color, eye_color) %>% mutate_at(vars(hair_color:eye_color), ~ gsub("none", NA, .)) %>% drop_na
starwars %>% summarise_all(~ sum(is.na(.)))
starwars %>% pmap_int(starwars, ~ sum(is.na(c(...))))
?pmap
pmap_int(starwars, ~ sum(is.na(c(...))))
keep <- pmap_int(starwars, ~ sum(is.na(c(...)))) < 4
sum(keep)
length(keep)
sw2 <- starwars %>% filter(keep)
sw2
dim(sw2)
as.matrix(starwars)
apply(starwars, 1, function(x) sum(is.na(x)))
apply(starwars, 2, function(x) sum(is.na(x)))
starwars %>% summarise_all(~ sum(is.na(.)))
is.na(starwars)
colSums(is.na(starwars))
rowSums(is.na(starwars))
