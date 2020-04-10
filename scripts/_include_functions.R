#!/usr/bin/env Rscript

# Equation for outputing linear model results
lm_eqn = function(df){
    m = lm(y ~ x, df)
    eq = substitute(italic(y) == b*italic(x) + a*","~~italic(r)^2~"="~r2, 
            list(   a = as.character(format(coef(m)[1], digits = 2)), 
                    b = as.character(format(coef(m)[2], digits = 2)), 
                    r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}
