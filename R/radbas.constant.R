"radbas.constant" <-
function (m, d) 
{
    if (d%%2 == 0) {
        Amd <- (((-1)^(1 + m + d/2)) * (2^(1 - 2 * m)) * (pi^(-d/2)))/(gamma(m) * 
            gamma(m - d/2 + 1))
    }
    else {
        Amd <- (gamma(d/2 - m) * (2^(-2 * m)) * (pi^(-d/2)))/gamma(m)
    }
    Amd
}
