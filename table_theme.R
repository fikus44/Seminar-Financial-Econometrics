
library("knitr")
library("kableExtra")


table_theme = function(data, colnames, caption) {
  
  data %>%
    kable(format = "latex", #Output er kompatibelt med LaTeX
          digits = 2,  # no. of digits after comma
          col.names = colnames, #søjlenavne
          align = c("l", rep("c", length(colnames) - 2), "r"), #Søjle alignement, første søjle er venstre aligned, herefter er alle søjle på nær den sidste søjle center aligned. Den sidste søjle er højre aligned
          caption = caption, #Tabel titel. Husk at caption skal angives som string i funktion. 
          format.args = list(big.mark = ",", scientific = FALSE), #Tilføjer  1.000 tals separator. Scientific = FALSE betyder vi IKKE benytter os af e til at denote tal
          escape = TRUE, #Output tager højde for special characters i LaTeX som f.eks. bliver _ til \_
          booktabs = TRUE, #Benytter os af booktabs package i LaTeX
          linesep =  '' #Intet ekstra line space hver x. linje. Kan ændres ved c("", "", "","","\\addlinespace")
            
          )
}


#Til input:

# Data skal være dataframe eller tibble
# colnames skal være vektor, hvor længden af vektor er lig med antallet af søjler i data
# caption skal være angivet som string

# Dokumentation:

# https://bookdown.org/yihui/rmarkdown-cookbook/kable.html
# https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf

# Tables kan customizes yderligere ved brug af bl.a. kable_styling() funktionen.
# Jeg vil foreslå, man ser på dokuemntation, hvis der er et specifikt customization
# problem som tabellen skal imødekomme. 

#Tester table_theme

data(mtcars)
test = head(mtcars)
names = c("mpg", "cyl", "disp", "hp", "drat", "wt", "qsec", "vs", "am", "gear", "carb")

table_theme(test, colnames = names, caption = "test af table_theme funktionen") %>%
  kable_styling(latex_options = "scale_down") #Scaler table til margins i LaTeX
