"WQSdi" <-
function (x, transpose = F) 
{
    if (!transpose) 
        t(WQSi(t(WQSi(x))))
    else {
        WQSi.T(t(WQSi.T(t(x))))
    }
}
