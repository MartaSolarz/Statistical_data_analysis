# Marta Solarz, Karolina Zdunek
# ZADANIE
# Analiza stanu drzew - klon i lipa
# Dane hiperspektralne

# Całość analizy wykonać dla porównania stanu (gorszy - bliżej drogi/lepszy - dalej od drogi) dla klonu i lipy oddzielnie
# Oraz w drugim wariancie wszystkie gatunki i stany łącznie (całość drzew, bez podziału na klasy)
# Czyli analiza dla 2 grup: (2 grupy: klon / lipa) i dla 4 grup (klon G/L, lipca G/L)
# 1. różnica istotnie statystystyczne między grupami w wartościach wskaźników wyliczonych dla tych drzew
#    wylicz 5 wskaźników teledetekcyjnych i czy się różnicują
# 2. przeprowadź korelację dla każdego gatunku oddzielnie (bez podziału na stan) dla kanałów
#    czyli wybierzemy kanały, które są najbardziej skorelowane i te najmniej
#    (nas w dalszym kroku interesują kanały najmniej skorelowane ze sobą)
# 3. wybierz z interesujących kanałów max. 10 (różne zakresy spektralne)
#    i sprawdź czy istnieją istotne statystycznie różnice pomiędzy badanymi gatunkami (wariant 4 typy)

# Każdy wynik interpretacji potwierdź graficznie, podaj interpretacje wyników, podaj hipotezy itp.

library(ggplot2)
library(dplyr)
library(reshape2)
library(dunn.test)

choose_numeric_columns <- function(df) {
  return(select_if(df, is.numeric))
}

find_no_gauss_variables <- function(df, alpha) { # return H1 variables
  dane <- select_if(df, function(x) {
    result <- shapiro.test(x)
    return(result$p.value < alpha)
  })
  return(names(dane))
}

find_variables_with_difference_median_two_groups <- function(df1, df2, alpha) { # return H1 variables
  no_diff <- c()
  for (i in seq(ncol(df1))) {
    result <- wilcox.test(df1[,i], df2[,i])
    if (result$p.value < alpha) {
      no_diff <- c(no_diff, colnames(df1)[i])
    }
  }
  return(no_diff)
}

find_variables_with_difference_median_four_groups <- function(num_cols, group_cols, alpha = 0.05) {  # return H1 variables
  no_diff <- c()
  for (i in seq(ncol(num_cols))) {
    result <- kruskal.test(num_cols[,i] ~ group_cols)
    if (result$p.value < alpha) {
      no_diff <- c(no_diff, colnames(num_cols)[i])
    }
  }
  return(no_diff)
}

setwd("/home/websensa/Studies/Programming_and_modeling_of_spatial_data/projekt")
dane <- read.table("dane_klon_lipa.csv", sep=";", dec = ",", header = TRUE)

# Sprawdzenie poprawności typów danych
str(dane)

# dodanie kolumny z połączonym gatunkiem i stanem
dane$rodzaj_stan <- paste(dane$rodzaj, dane$stan, sep = "_")

# --------------------------------------------------------------
# ZADANIE 1
# --------------------------------------------------------------
# Wybrane wskaźniki NDVI, WBI, NDII, MSI, NDVI705 pozwalające na ocenę stanu kondycji roślin
# Formuły:
ndvi <- as.numeric((dane$B121 - dane$B83)/(dane$B121 + dane$B83))
wbi <- as.numeric(dane$B167/dane$B153)
ndii <- as.numeric((dane$B122 - dane$B283)/(dane$B122 + dane$B283))
msi <- as.numeric(dane$B283/dane$B122)
ndvi705 <- as.numeric((dane$B106 - dane$B90)/(dane$B106 + dane$B90))

df <- data.frame(
  Gatunek = dane$rodzaj,
  Stan = dane$stan,
  Rodzaj_stan = dane$rodzaj_stan,
  NDVI = ndvi,
  WBI = wbi,
  NDII = ndii,
  MSI = msi,
  NDVI705 = ndvi705
)

lipa <- choose_numeric_columns(df[df$Gatunek == "lipa",])
klon <- choose_numeric_columns(df[df$Gatunek == "klon",])
lipa_L <- choose_numeric_columns(subset(df, df$Rodzaj_stan == "lipa_L"))
lipa_G <- choose_numeric_columns(subset(df, df$Rodzaj_stan == "lipa_G"))
klon_L <- choose_numeric_columns(subset(df, df$Rodzaj_stan == "klon_L"))
klon_G <- choose_numeric_columns(subset(df, df$Rodzaj_stan == "klon_G"))

# --------------------------------------------------------------
# WIZUALIZACJE
# Wykresy wskaźników w zależności od gatunku i stanu
ggplot(data = df, aes(x = Stan, y = NDVI, color = Gatunek)) +
  geom_boxplot()
ggplot(data = df, aes(x = Stan, y = NDII, color = Gatunek)) +
  geom_boxplot()
ggplot(data = df, aes(x = Stan, y = NDVI705, color = Gatunek)) +
  geom_boxplot()
ggplot(data = df, aes(x = Stan, y = WBI, color = Gatunek)) +
  geom_boxplot()
ggplot(data = df, aes(x = Stan, y = MSI, color = Gatunek)) +
  geom_boxplot()

# Wykresy gęstości
# dla lipy
par(mfrow=c(3,2))
for (i in seq(ncol(lipa))) {
  par(mar = c(2,2,1,1))
  density_lipa <- density(lipa[,i], na.rm=TRUE)
  density_lipa_L <- density(lipa_L[,i], na.rm=TRUE)
  density_lipa_G <- density(lipa_G[,i], na.rm=TRUE)
  max_y_lim <- max(density_lipa$y, density_lipa_L$y, density_lipa_G$y)
  plot(density_lipa, col = "black", lwd = 1, lty = 1, ylim=c(0, max_y_lim), main="", xaxt = "n")
  mtext(side = 1, text= names(lipa)[i], line=1)
  lines(density_lipa_L, col = "green", lwd = 1, lty = 2)
  lines(density_lipa_G, col = "red", lwd = 1, lty = 2)
}
plot.new()
legend("center", legend = c("Lipa", "Lipa lepsza", "Lipa gorsza"),
         col = c("black", "green", "red"), lwd = 1, lty = c(1,2,2))

# dla klonu
par(mfrow=c(3,2))
for (i in seq(ncol(lipa))) {
  par(mar = c(2,2,1,1))
  density_klon <- density(klon[,i], na.rm=TRUE)
  density_klon_L <- density(klon_L[,i], na.rm=TRUE)
  density_klon_G <- density(klon_G[,i], na.rm=TRUE)
  max_y_lim <- max(density_klon$y, density_klon_L$y, density_klon_G$y)
  plot(density_klon, col = "black", lwd = 1, lty = 1, ylim=c(0, max_y_lim), main="", xaxt = "n")
  mtext(side = 1, text= names(lipa)[i], line=1)
  lines(density_klon_L, col = "green", lwd = 1, lty = 2)
  lines(density_klon_G, col = "red", lwd = 1, lty = 2)
}
plot.new()
legend("center", legend = c("Klon", "Klon lepszy", "Klon gorszy"),
         col = c("black", "green", "red"), lwd = 1, lty = c(1,2,2))

# --------------------------------------------------------------
# TESTOWANIE
# Test normalności
# H0: rozkład wskaźnika X dla grupy Y jest rozkładem Gaussa
# H1: rozkład wskaźnika X dla grupy Y jest różny od rozkładu Gaussa
# Wybieramy poziom istotności alpha = 0.05
find_no_gauss_variables(lipa, 0.05) # --> brak wskaźników wykazujących rozkład normalny, dla wszystkich przyjmujemy H1
find_no_gauss_variables(klon, 0.05) # --> brak wskaźników wykazujących rozkład normalny, dla wszystkich przyjmujemy H1
find_no_gauss_variables(lipa_L, 0.05) # --> brak wskaźników wykazujących rozkład normalny, dla wszystkich przyjmujemy H1
find_no_gauss_variables(lipa_G, 0.05) # --> brak wskaźników wykazujących rozkład normalny, dla wszystkich przyjmujemy H1
find_no_gauss_variables(klon_L, 0.05) # --> brak wskaźników wykazujących rozkład normalny, dla wszystkich przyjmujemy H1
find_no_gauss_variables(klon_G, 0.05) # --> w przypadku "NDII" oraz "MSI" brak podstaw do odrzucenia H0, dla pozostałych przyjmujemy H1

# Wnioski:
# A zatem dla każdego wskaźnika z każdej z analizowanych grup
# (1. lipa vs klon; 2. lipa lepsza vs lipa gorsza; 3. klon lepszy vs klon gorszy)
# conajmniej jedna z grup ma rozkład danych statystycznie istotnie różny od normalnego
# zatem dalej musimy stosować do wszystkich grup testy nieparametryczne.

# --------------------------------------------------------------
# Test istotnej statystycznie różnicy median między grupami
# poziom istotności alpha = 0.05
# H0: brak istotnych statystycznie różnic dla wskaźnika X dla grup A i B
# H1: istnieją statystczynie istotne różnice dla wskaźnika X dla grup A i B

find_variables_with_difference_median_two_groups(lipa, klon, 0.05) # --> dla wszystkich wskaźników przyjmujemy H1
find_variables_with_difference_median_two_groups(lipa_L, lipa_G, 0.05) # --> brak podstaw do odrzucenia H0 dla MSI, WBI i NDII, dla NDVI i NDVI705 przyjmujemy H1
find_variables_with_difference_median_two_groups(klon_L, klon_G, 0.05) # --> brak podstaw do odrzucenia H0 dla MSI, WBI i NDII, dla NDVI i NDVI705 przyjmujemy H1

# Wnioski:
# Zauważalne różnice między wszystkimi wskaźnikami przy porównaniu lipa vs klon
# Dla obu gatunków widoczne są statystycznie istotne różnice w wartości NDVI i NDVI705 przy porównaniu stanu gorszego vs lepszego

# Test dla 4 grup razem
find_variables_with_difference_median_four_groups(as.data.frame(choose_numeric_columns(df)), df$Rodzaj_stan) # --> dla wszystkich wskaźników przyjmujemy H1

# ponieważ dla wszystkich przypadków otrzymaliśmy p < alpha należy przeprowadzić test post-hoc
# wybieramy metodę Holma ze względu na średni poziom konserwatywności metody
dunn.test(df$NDVI, df$Rodzaj_stan, method = "holm") # istotna statystycznie różnca między każdą z grup
dunn.test(df$NDII, df$Rodzaj_stan, method = "holm") # istotna statystycznie różnica między klon_G i lipa_G, lipa_G i klon_L, lipa_L i klon_G, klon_L i lipa_L
dunn.test(df$MSI, df$Rodzaj_stan, method = "holm") # istotna statystycznie różnica między klon_G i lipa_G, lipa_G i klon_L, lipa_L i klon_G, klon_L i lipa_L
dunn.test(df$WBI, df$Rodzaj_stan, method = "holm") # istotna statystycznie różnica między lipa_G i klon_L, klon_L i lipa_L
dunn.test(df$NDVI705, df$Rodzaj_stan, method = "holm") # istotna statystycznie różnca między każdą z grup

# Wnioski:
# Zauważalne różnice między wszystkimi wskaźnikami przy porównaniu lipa vs klon
# Dla porównania 4 grup widoczne są istotne statystycznie różnice dla wszystkich wskaźników

# --------------------------------------------------------------
# ZADANIE 2
# --------------------------------------------------------------
# Przeprowadź korelację dla każdego gatunku oddzielnie (bez podziału na stan) dla kanałów
#    czyli wybierzemy kanały, które są najbardziej skorelowane i te najmniej
#    (nas w dalszym kroku interesują kanały najmniej skorelowane ze sobą)

# wybór kolumn do analizy
wybrane_kolumny <- data.frame(dane[9:ncol(dane)])
wybrane_kolumny[, 2:431] <- apply(wybrane_kolumny[, 2:431], 2, as.numeric)
str(wybrane_kolumny)

# wybranie tylko tych co mają charkter numeryczny w podziale na rodzaj
kanaly_lipa <- choose_numeric_columns(wybrane_kolumny[wybrane_kolumny$rodzaj == "lipa",])
kanaly_klon <- choose_numeric_columns(wybrane_kolumny[wybrane_kolumny$rodzaj == "klon",])

# TESTOWANIE KORELACJI
# Aby przejść do testowania korelacji między kanałami dla gatunków musimy wybrać test
# w zależności od rozkładu wartości w poszczególnych kanałach
# Testowanie normalności
# H0: rozkład wartości kanału X dla grupy Y jest rozkładem Gaussa
# H1: rozkład wartości kanału X dla grupy Y jest różny od rozkładu Gaussa
# Wybieramy poziom istotności alpha = 0.05

kanaly_lipa_H1 <- find_no_gauss_variables(kanaly_lipa, 0.05)
kanaly_klon_H1 <- find_no_gauss_variables(kanaly_klon, 0.05)

# wybieramy zmienne o rozkładzie normalnym
kanaly_lipa_H0 <- setdiff(names(kanaly_lipa), kanaly_lipa_H1)
kanaly_klon_H0 <- setdiff(names(kanaly_klon), kanaly_klon_H1)

df_lipa_gauss <- kanaly_lipa[kanaly_lipa_H0]
df_klon_gauss <- kanaly_klon[kanaly_klon_H0]

df_lipa_nie_gauss <- kanaly_lipa[kanaly_lipa_H1]
df_klon_nie_gauss <- kanaly_klon[kanaly_klon_H1]

# Korelacje
# Funkcje pomocnicze:
get_new_p_matrix <- function (x) {
  return(matrix(nrow=ncol(x), ncol=ncol(x)))
}

get_significant <- function(matrix_p, matrix_corr, alpha) {
  significant <- matrix_p
  significant[matrix_p < alpha] <- 1
  significant[matrix_p >= alpha] <- NA

  significant_matrix <- matrix_corr
  significant_matrix[is.na(significant)] <- NA

  return(significant_matrix)
}

get_corr_p_value <- function(x, matrix_p, matrix_corr, corr_method, alpha) {
  for (i in seq(ncol(x))) {
    for (j in seq(ncol(x))) {
      cor_val <- cor.test(x[, i], x[, j],
                                use="complete.obs", method = corr_method, exact = FALSE)
      p_val <- cor_val$p.value
      matrix_p[i,j] <- p_val
    }
  }
  return(get_significant(matrix_p, matrix_corr, alpha))
}

create_graphics <- function(matrix_corr, title) {
  ggplot(data = reshape2::melt(matrix_corr), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red") +
  # geom_text(aes(label = round(value, 1))) +  # można odkomentować jeśli chcemy mieć etykiety
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Kanały", y = "Kanały", fill = "Wsp. korelacji", title = title)
}

# H0: brak zależności między kanałami X i Y
# H1: istnieje zależność między kanałami X i Y
# alpha = 0.05 - domyślny

# LIPA
# a) kanały o rozkładzie normalnym
macierz_korelacji_lipa_gauss <- cor(df_lipa_gauss, use = "pairwise.complete.obs", method = "pearson")
macierz_p_lipa_gauss <- get_new_p_matrix(df_lipa_gauss)
macierz_istotnych_korelacji_lipa_gauss <- get_corr_p_value(df_lipa_gauss, macierz_p_lipa_gauss,
  macierz_korelacji_lipa_gauss, "pearson", 0.05)
create_graphics(macierz_istotnych_korelacji_lipa_gauss,
                "Istotna statystycznie korelacja dla lipy - kanały o rozkładzie normalnym")

# Wnioski:
# kanały mocno skorelowane ze sobą to B1:B4 - B1:B4, B89-B89, B152:B163 (bez B159, B160) - B152:B163 (bez B159, B160), B164:B248 - B164:B248,
# kanały słabo skorelowane to B1:B4 - B152:B248 (najniższa korelacaja)
# Ogólna zależność -> kanały położone blisko siebie są bardziej skorelowane, niż te które dzieli większa odległość, duża korelacja obserwowana po przekątnej

# b) kanały o rozkładzie różnym od normalnego
macierz_korelacji_lipa_nie_gauss <- cor(df_lipa_nie_gauss, use = "pairwise.complete.obs", method = "spearman")
macierz_p_lipa_nie_gauss <- get_new_p_matrix(df_lipa_nie_gauss)
macierz_istotnych_korelacji_lipa_nie_gauss <- get_corr_p_value(df_lipa_nie_gauss, macierz_p_lipa_nie_gauss,
  macierz_korelacji_lipa_nie_gauss, "spearman", 0.05)
create_graphics(macierz_istotnych_korelacji_lipa_nie_gauss,
                "Istotna statystycznie korelacja dla lipy - kanały o rozkładzie różnym od normalnego")

# Wnioski:
# kanały mocno skorelowane ze sobą to B5:B88 - B5:B88,  B90:B160 (bez B152-B158) -  B90:B160 (bez B152-B158), B166 - B197 (bez B168-B169), B249:430 - B249:430
# kanały słabo skorelowane to  B5:B88 - B166 - B197 (bez B168-B169)
# Ogólna zależność -> kanały położone blisko siebie są bardziej skorelowane, niż te które dzieli większa odległość, duża korelacja obserwowana po przekątnej

# KLON
# a) kanały o rozkładzie normalnym
macierz_korelacji_klon_gauss <- cor(df_klon_gauss, use = "pairwise.complete.obs", method = "pearson")
macierz_p_klon_gauss <- get_new_p_matrix(df_klon_gauss)
macierz_istotnych_korelacji_klon_gauss <- get_corr_p_value(df_klon_gauss, macierz_p_klon_gauss,
  macierz_korelacji_klon_gauss, "pearson", 0.05)
create_graphics(macierz_istotnych_korelacji_klon_gauss,
                "Istotna statystycznie korelacja dla klonu - kanały o rozkładzie normalnym")

# Wnioski:
# kanały mocno skorelowane ze sobą to B36:B52 - B36:B52, kanały z zakresu B114:B164 - B114:B164, B317:B326 - B317:B326,
# kanały słabo skorelowane to B36:B52 - B114:B164
# Ogólna zależność -> kanały położone blisko siebie są bardziej skorelowane, niż te które dzieli większa odległość, duża korelacja obserwowana po przekątnej


# b) kanały o rozkładzie różnym od normalnego
macierz_korelacji_klon_nie_gauss <- cor(df_klon_nie_gauss, use = "pairwise.complete.obs", method = "spearman")
macierz_p_klon_nie_gauss <- get_new_p_matrix(df_klon_nie_gauss)
macierz_istotnych_korelacji_klon_nie_gauss <- get_corr_p_value(df_klon_nie_gauss, macierz_p_klon_nie_gauss,
  macierz_korelacji_klon_nie_gauss, "spearman", 0.05)
create_graphics(macierz_istotnych_korelacji_klon_nie_gauss,
                "Istotna statystycznie korelacja dla klonu - kanały o rozkładzie różnym od normalnego")

# Wnioski:
# kanały mocno skorelowane ze sobą to B1:B35 - B1:B35, B35:B88 - B35:B88, B89:B113 - B89:B113, B164:B430 - B164:B430
# kanały słabo skorelowane to B1:B35 - B166:B197
# Ogólna zależność -> kanały położone blisko siebie są bardziej skorelowane, niż te które dzieli większa odległość,  duża korelacja obserwowana po przekątnej

# UWAGA: aby zwiększyć czytelność wykresów można generować wykresy dla części kanałów np:
# create_graphics(macierz_korelacji_lipa_nie_gauss[1:100, 1:100], "Korelacja dla 100 kanałów")

# --------------------------------------------------------------
# ZADANIE 3
# --------------------------------------------------------------
# Wybierz z interesujących kanałów (jak najmniejsza korelacja) max. 10 (różne zakresy spektralne)
# i sprawdź czy istnieją istotne statystycznie różnice pomiędzy badanymi gatunkami (wariant 4 typy)

# Wybieramy mało skorelowane kanały z różnych kanałów spektralnych do dalszej analizy
kanaly_do_analizy <- dane[c("rodzaj", "stan", "rodzaj_stan", "B4", "B50", "B61", "B117", "B164", "B217", "B240", "B326", "B353")]

lipa_l <- kanaly_do_analizy[kanaly_do_analizy$rodzaj_stan == "lipa_L",]
lipa_g <- kanaly_do_analizy[kanaly_do_analizy$rodzaj_stan == "lipa_G",]
klon_l <- kanaly_do_analizy[kanaly_do_analizy$rodzaj_stan == "klon_L",]
klon_g <- kanaly_do_analizy[kanaly_do_analizy$rodzaj_stan == "klon_G",]

# Testowanie normalności
# H0: rozkład wartości kanału X dla grupy Y jest rozkładem Gaussa
# H1: rozkład wartości kanału X dla grupy Y jest różny od rozkładu Gaussa
# Wybieramy poziom istotności alpha = 0.05
find_no_gauss_variables(choose_numeric_columns(lipa_l), 0.05) # "B4"   "B50"  "B61"  "B117" "B164" "B217" "B240" "B326" "B353"
find_no_gauss_variables(choose_numeric_columns(lipa_g), 0.05) # "B50"         "B61"                              "B326" "B353"
find_no_gauss_variables(choose_numeric_columns(klon_l), 0.05) # "B4"   "B50"         "B117" "B164" "B217" "B240" "B326" "B353"
find_no_gauss_variables(choose_numeric_columns(klon_g), 0.05) # "B4"          "B61"                "B217" "B240"        "B353"

# A zatem zauważmy, że zarówno dla klonu, jak i dla lipy i dla wszystkich wybranych kanałów
# przynajmniej jeden kanał z pary x_gorsza i x_lepsza ma rozkład różny od normalnego
# zatem korzystamy z testów nieparametrycznych

# Test istotnej statystycznie różnicy median między grupami
# poziom istotności alpha = 0.05
# H0: brak istotnych statystycznie różnic dla wskaźnika X dla grup A i B
# H1: istnieją statystczynie istotne różnice dla wskaźnika X dla grup A i B

# między 2 grupami
find_variables_with_difference_median_two_groups(choose_numeric_columns(lipa_l), choose_numeric_columns(lipa_g), 0.05) #        "B50"  "B61"  "B117" "B164" "B217" "B240" "B326"
find_variables_with_difference_median_two_groups(choose_numeric_columns(klon_l), choose_numeric_columns(klon_g), 0.05) # "B4"   "B50"  "B61"  "B117" "B164" "B217" "B240" "B326"

# A zatem dla lipy istnieje istotna statystycznie różnica między stanami dla kanałów:
# "B50"  "B61"  "B117" "B164" "B217" "B240" "B326"
# a dla kanałów "B4" "B353" nie jest ona zauważalna.

# A dla klonu istnieje istotna statystycznie różnica między stanami dla kanałów:
# "B4"   "B50"  "B61"  "B117" "B164" "B217" "B240" "B326"
# a dla kanału "B353" nie jest ona zauważalna.

# miedzy 4 grupami
find_variables_with_difference_median_four_groups(as.data.frame(choose_numeric_columns(kanaly_do_analizy)), kanaly_do_analizy$rodzaj_stan) # "B4"   "B50"  "B61"  "B117" "B164" "B217" "B240" "B326" "B353"
# A zatem istnieje istotna statystycznie różnica między stanami dla wszystkich wybranych kanałów
# Test post-hoc:
dunn.test(kanaly_do_analizy$B4, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między    klon_G i klon_L, klon_G i lipa_G, klon_G i lipa_L, klon_L i lipa_G
dunn.test(kanaly_do_analizy$B50, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między   klon_G i klon_L,                  klon_G i lipa_L, klon_L i lipa_G,                  lipa_G i lipa_L
dunn.test(kanaly_do_analizy$B61, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica dla wszystkich kanałów
dunn.test(kanaly_do_analizy$B117, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między  klon_G i klon_L,                  klon_G i lipa_L, klon_L i lipa_G, klon_L i lipa_L, lipa_G i lipa_L
dunn.test(kanaly_do_analizy$B164, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między  klon_G i klon_L,                  klon_G i lipa_L, klon_L i lipa_G, klon_L i lipa_L, lipa_G i lipa_L
dunn.test(kanaly_do_analizy$B217, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między  klon_G i klon_L, klon_G i lipa_G,                  klon_L i lipa_G, klon_L i lipa_L, lipa_G i lipa_L
dunn.test(kanaly_do_analizy$B240, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między  klon_G i klon_L, klon_G i lipa_G,                  klon_L i lipa_G, klon_L i lipa_L, lipa_G i lipa_L
dunn.test(kanaly_do_analizy$B326, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica dla wszystkich kanałów
dunn.test(kanaly_do_analizy$B353, kanaly_do_analizy$rodzaj_stan, method = "holm") # istotna statystycznie różnica między                   klon_G i lipa_G, klon_G i lipa_L, klon_L i lipa_G, klon_L i lipa_L

num <- choose_numeric_columns(kanaly_do_analizy)
par(mfrow=c(3,3))
for (i in seq_along(num)) {
  par(mar = c(2,2,1,1))
  density <- density(num[,i], na.rm=TRUE)
  density_lipa <- density(subset(num[,i], kanaly_do_analizy$rodzaj == 'lipa'), na.rm=TRUE)
  density_klon <- density(subset(num[,i], kanaly_do_analizy$rodzaj == 'klon'), na.rm=TRUE)
  max_y_lim <- max(density_lipa$y, density$y, density_klon$y)
  plot(density, col = "black", lwd = 1, lty = 1, ylim=c(0, max_y_lim), main="", xaxt = "n")
  mtext(side = 1, text= colnames(num)[i], line=1)
  lines(density_lipa, col = "green", lwd = 1, lty = 2)
  lines(density_klon, col = "red", lwd = 1, lty = 2)
  legend("right", legend = c("klon+lipa", "lipa", "klon"),
         col = c("black", "green", "red"), lwd = 1, lty = c(1,2,2))
}

# Boxploty
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B4, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B4 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B50, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B50 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B61, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B61 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B117, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B117 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B164, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B164 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B217, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B217 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B240, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B240 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B326, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B326 w podziale na badane gatunki")
ggplot(data = kanaly_do_analizy, aes(x = rodzaj_stan, y = B353, color = rodzaj_stan)) +
  geom_boxplot() + ggtitle("Rozkład wartości kanału B353 w podziale na badane gatunki")
